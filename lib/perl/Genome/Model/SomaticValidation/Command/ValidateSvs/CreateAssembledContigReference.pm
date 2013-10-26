package Genome::Model::SomaticValidation::Command::ValidateSvs::CreateAssembledContigReference;

use strict;
use warnings;

use Genome;

use POSIX qw(ceil);

class Genome::Model::SomaticValidation::Command::ValidateSvs::CreateAssembledContigReference {
    is => 'Command::V2',
    has_input => [
        build_id => {
            is => 'Text',
            doc => 'ID of the somatic validation build',
            is_output => 1,
        },
        skip => {
            is => 'Boolean',
            doc => 'signal from first step whether or not to run',
            default_value => 0,
        },
    ],
    has => [
        build => {
            is => 'Genome::Model::Build::SomaticValidation',
            id_by => 'build_id',
        },
        output_dir => {
            is_output => 1,
            is => 'Text',
            doc => 'Place where the output goes',
            is_calculated => 1,
            calculate_from => 'build',
            calculate => q{ return join("/", $build->data_directory, "/validation/sv"); }
        },
    ],
    has_output => [
        reference_build_id => {
            is => 'Text',
            doc => 'ID of the newly create reference sequence',
            is_optional => 1,
        },
    ],
    doc => 'creates a reference sequence for validating by alignment',
};

sub sub_command_category { 'pipeline steps' }

sub execute {
    my $self = shift;
    my $build = $self->build;

    if($self->skip) {
        $self->status_message("skip signal received. not running.");
        return 1;
    }

    #input parameters
    my $calls = join("/", $self->output_dir, "assembly_output.csv.merged");
    my $sample_id = $build->tumor_sample->name;

    #put calls into a hash for use when parsing fasta input file
    my %calls;
    my $calls_fh = Genome::Sys->open_file_for_reading($calls);
    while (my $line = $calls_fh->getline) {
        if ($line =~ /^#/) {
            #check header
            unless ($line =~ /^#ID\tCHR1\tOUTER_START\tINNER_START\tCHR2/) {
                die $self->error_message("SV calls do not seem to be in merged format with header #ID\tCHR1\tOUTER_START\tINNER_START\tCHR2...");
            }
            next;
        }

        my ($id) = split("\t",$line);
        $calls{$id}++;
    }
    $calls_fh->close;

    my $contigs_file = $self->_process_fasta(\%calls);

    #get refseq build id from model unless already defined
    my $ref_seq_build = $build->reference_sequence_build;

    #create reference sequence using new contigs (define new reference and track new reference build)
    my $ref_model_id = $ref_seq_build->model_id;

    my $version = "500bp_assembled_contigs_sv";
    #don't overwrite an existing model...
    my $prefix = $sample_id . "_SV_Contigs";
    $version = $self->check_ref_build_name($prefix, $version);

    my $new_ref_cmd = Genome::Model::Command::Define::ImportedReferenceSequence->create(
        species_name => 'human',
        use_default_sequence_uri => '1',
        derived_from => $ref_seq_build,
        append_to => $ref_seq_build,
        version => $version,
        fasta_file => $contigs_file,
        prefix => $prefix,
        server_dispatch => 'inline',
        is_rederivable => 1,
    );
    unless ($new_ref_cmd->execute) {
        $self->error_message('Failed to execute the definition of the new reference sequence with added contigs.');
        return;
    }
    my $new_ref_build_id = $new_ref_cmd->result_build_id;
    my $new_ref_build = Genome::Model::Build->get($new_ref_build_id);
    unless ($new_ref_build->status eq 'Succeeded') {
        die $self->error_message('New reference build not successful.');
    }

    $self->reference_build_id($new_ref_build_id);
    return 1;
}

sub _process_fasta {
    my $self = shift;
    my $calls = shift;

    my $fasta = join("/", $self->output_dir, "assembly_output.fasta.merged");
    my $contigs_file = join("/", $self->output_dir, "contigs");
    my $review_file = join("/", $self->output_dir, "manual_review");

    my $fasta_fh = Genome::Sys->open_file_for_reading($fasta);
    my $contig_fh = Genome::Sys->open_file_for_writing($contigs_file);
    my $review_fh = Genome::Sys->open_file_for_writing($review_file);

    #print header for manual review file
    print $review_fh "#ID\tChr1\tPos1\tChr2\tPos2\tType\tContig_Length\tBreakpoint_Estimation\n";

    #parse and relabel fasta contigs associated with the calls in the sv_calls_file
    #sample line:
    #>1.10,LUC12_SV,Var:1.93179141.1.93179176.DEL.33.+-,Ins:336-338,Length:631,KmerCoverage:45.32,Strand:+,Assembly_Score:228.96,PercNonRefKmerUtil:9,Ref_start:93178807,Ref_end:93179470,Contig_start:1,Contig_end:631,TIGRA
    #ccttgcatatttttaagttgacatctacaatttttcaccataagtttaaatagttgcaaa

    {
        my ($id,$source,$var,$ins,$length,$kmer,$strand,$score,$nonrefkmer,$ref_start,$ref_end);
        LINE: while (my $line = $fasta_fh->getline) {
            chomp $line;

            if ($line =~ /^\>/) {
                ($id,$source,$var,$ins,$length,$kmer,$strand,$score,$nonrefkmer,$ref_start,$ref_end) = split ",",$line;

                #make sure id is from a somatic site
                $id =~ s/^>//;
                unless ($calls->{$id}) {
                    undef($id);
                    next LINE;
                }
            }

            if(not $id) {
                #part of a contig we don't care about
                next LINE;
            }

            #clean up some variables from the assembly header
            $var =~ s/^Var:([\w\d\.]+)\.\d+\.[+-]+$/$1/;
            (undef,$length) = split(":",$length);
            (undef,$strand) = split(":",$strand);
            (undef,$ins) = split(":",$ins);

            #find breakpoint within contig sequence
            my ($microhomology_start, $microhomology_stop) = $ins =~ /(\d+)\-(\d*)/;
            my $breakpoint = $microhomology_start;
            #stop position is not always given. if it is, split the difference (should be a small difference)
            if ($microhomology_stop) { $breakpoint = ceil(($microhomology_start + $microhomology_stop)/2); }

            #grab more lines until contig sequence is fully obtained
            my $contig = '';
            my $contig_line;
            while ($contig_line = $fasta_fh->getline and $contig_line !~ /^>/) {
                chomp $contig_line;
                $contig .= $contig_line;
            }

            unless($contig) {
                die ($self->error_message("Contig had no lines"));
            }

            #handle negative stranded contigs
            if ($strand eq '-') {
                $contig =~ tr/ACGTacgt/TGCAtgca/;
                $contig = reverse($contig);
                $breakpoint = $length - $breakpoint;
            }

            #print output contig
            my $description = join(" ","ID_".$id."_CALL_".$var,"Breakpoint_Within_Contig:~".$breakpoint,"Contig_Length:".$length,"Original_".$score,"Original_Contig_Strand:".$strand) . "\n";
            print $contig_fh ">",$description;
            while ($contig) { print $contig_fh substr($contig,0,80,""),"\n"; }

            #print output manual review file
            $var =~ s/\./\t/g;
            my $review_line = join("\t",$id,$var,$length,$breakpoint);
            print $review_fh "$review_line\n";

            if($contig_line and $contig_line =~ /^>/) {
                #we ran into the next contig... process it now.
                $line = $contig_line;
                redo LINE;
            }
        }
    }

    $fasta_fh->close;
    $contig_fh->close;
    $review_fh->close;

    return $contigs_file;
}

sub check_ref_build_name {
    my $self = shift;
    my $sample_id = shift;
    my $version = shift;

    #FIXME this should be all about a specific build type...
    my @builds = Genome::Model::Build->get("model.name LIKE" => $sample_id . "-human%");

    my $max = -1;
    for my $build (@builds) {
        my $v = $build->version;
        #if we have models with a suffix already, store the highest suffix
        if ($v =~ /$version-(\d+)/){
            if ($1 > $max){
                $max = $1;
            }
            #else if we have a match at all for this model name
        } 
        elsif ($v =~ /$version/) {
            my $ver = 0;
            $max = $ver if $ver > $max;
        }
    }
    if ($max > -1) {
        $max++;
        $version = $version . "-$max";
    }
    return $version;
}

1;
