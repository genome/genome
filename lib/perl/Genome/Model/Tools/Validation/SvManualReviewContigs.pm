package Genome::Model::Tools::Validation::SvManualReviewContigs;

use strict;
use warnings;

use IO::File;
use Genome;
use POSIX;

class Genome::Model::Tools::Validation::SvManualReviewContigs {
    is => 'Command',
    has_input => [
        merged_sv_calls_file => {
            is => 'String',
            doc => 'path to SV calls that should be reviewed using contig alignments (in "merged" format)',
        },
        merged_assembly_fasta_file => {
            is => 'String',
            doc => 'path to assembly output fasta file associated with sv_calls_file (in "merged" format)',
        },
        contigs_output_file => {
            is => 'String',
            doc => 'output file of contigs which will be appended to reference sequence for realignment',
        },
        manual_review_output_file => {
            is => 'String',
            doc => 'contig alignment summary file to be used for manual review',
        },
        sample_identifier => {
            is => 'String',
            doc => 'some string to use in model names, etc, such as "BRC2"',
        },
    ],
    has_optional_input => [
        somatic_validation_build_id => {
            is => 'Number',
            doc => 'somatic-validation build ID (contains both tumor and normal)',
        },
        somatic_validation_model_id => {
            is => 'Number',
            is_deprecated => 1,
            doc => 'somatic-validation model ID (contains both tumor and normal). will grab last succeeded build and do some possibly inconsistent analysis with it',
        },
        tumor_val_model_id => {
            is => 'Number',
            doc => 'tumor validation model to copy for alignment to new reference', #FIXME This is also possibly bad--use builds.
        },
        normal_val_model_id => {
            is => 'Text',
            doc => 'normal validation model to copy for alignment to new reference', #FIXME This is also possibly bad--use builds.
        },
        reference_sequence_build_id => {
            is => 'Text',
            doc => 'reference sequence build to be used for creating new reference with appended contigs',
        },
    ],
    has_transient_optional_output => [
        realignment_model_ids => {
            is => 'Text',
            is_many => 1,
            doc => 'the ids of the models that we made for aligning to our new contigs',
        },
    ],
    doc => 'create and align validation data to an SV contig reference',
};

sub help_detail {
    return <<EOS
This tool takes as input a file of SV-calls and an associated assembly fasta output file. Assembly contigs are matched to each event, and then appendended to the reference sequence specified via --reference-sequence-build-id. This modified reference sequence is imported into the system, and then current validation builds are copied, with the copies modified from the original to align to the new reference sequence. A command is given via STDOUT for use in running builds of these realignment models, after which alignments to the SV contigs may be used for manual review.
EOS
}

sub execute {

    my $self = shift;

    #parse input parameters
    my $calls = $self->merged_sv_calls_file;
    my $fasta = $self->merged_assembly_fasta_file;
    my $contigs_file = $self->contigs_output_file;
    my $review_file = $self->manual_review_output_file;
    my $sample_id = $self->sample_identifier;

    my $input_model_type;

    if(defined($self->somatic_validation_build_id or $self->somatic_validation_model_id)){
        $input_model_type = "somval";
    } elsif(defined($self->tumor_val_model_id) && defined($self->normal_val_model_id)){
        $input_model_type = "pairedref";
    } else {
        $self->error_message("you mush specify either a somatic_validation_model_id or tumor_val_model_id, normal_val_model_id, and reference_sequence_build_id");
        return 0;
    }

    #put calls into a hash for use when parsing fasta input file
    my %calls;
    my $calls_fh = new IO::File $calls,"r";
    while (my $line = $calls_fh->getline) {

        #check header
        if ($line =~ /^#/) {
            unless ($line =~ /^#ID\tCHR1\tOUTER_START\tINNER_START\tCHR2/) {
                $self->error_message("SV calls do not seem to be in merged format with header #ID\tCHR1\tOUTER_START\tINNER_START\tCHR2...");
                return;
            }
            next;
        }

        #sample line:
        #ID     CHR1    OUTER_START     INNER_START     CHR2    INNER_END       OUTER_END       TYPE    ORIENTATION     MINSIZE MAXSIZE SOURCE  SCORES  Copy_Number
        #1.3     1       18739525        18739525        1       18739874        18739874        INV     ++      350     350     LUC12_SV        237

        chomp $line;
        my ($id,$chr1,undef,undef,$chr2) = split /\t/,$line;
        $calls{$id}++;

    }
    $calls_fh->close;

    #open fasta file and output files
    my $fasta_fh = new IO::File $fasta,"r";
    my $line = $fasta_fh->getline;
    my $contig_fh = new IO::File $contigs_file,"w";
    my $review_fh = new IO::File $review_file,"w";

    #print header for manual review file
    print $review_fh "#ID\tChr1\tPos1\tChr2\tPos2\tType\tContig_Length\tBreakpoint_Estimation\n";

    #check that ->eof will allow one to enter the loop as expected
    if ($fasta_fh->eof) { die "fasta file only had one line!!\n"; }

    #parse and relabel fasta contigs associated with the calls in the sv_calls_file
    while (! $fasta_fh->eof) {
        chomp $line;

        #sample line:
        #>1.10,LUC12_SV,Var:1.93179141.1.93179176.DEL.33.+-,Ins:336-338,Length:631,KmerCoverage:45.32,Strand:+,Assembly_Score:228.96,PercNonRefKmerUtil:9,Ref_start:93178807,Ref_end:93179470,Contig_start:1,Contig_end:631,TIGRA
        #ccttgcatatttttaagttgacatctacaatttttcaccataagtttaaatagttgcaaa

        if ($line =~ /^\>/) {
            my ($id,$source,$var,$ins,$length,$kmer,$strand,$score,$nonrefkmer,$ref_start,$ref_end) = split ",",$line;

            #make sure id is from a somatic site
            $id =~ s/^\>(.+)$/$1/;
            unless ($calls{$id}) { 
                $line = $fasta_fh->getline;
                next;
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

            #grab next line which should contain the first bit of sequence for this contig
            my $contig = $fasta_fh->getline;
            chomp $contig;

            #grab more lines until contig sequence is fully obtained
            $line = $fasta_fh->getline;
            while ($line !~ /^\>/) {
                chomp $line;
                $contig .= $line;
                if ($fasta_fh->eof) { last; } 
                else { $line = $fasta_fh->getline; }
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

        }

        else { $line = $fasta_fh->getline; }

    }
    $fasta_fh->close;
    $contig_fh->close;
    $review_fh->close;


    #-------------------------------------------------------
    #get refseq build id from model unless already defined
    my $ref_seq_build;
    if(defined($self->reference_sequence_build_id)){
        $ref_seq_build = Genome::Model::Build->get($self->reference_sequence_build_id);
    } else {
        if($input_model_type eq "somval"){
            my $build_or_model = Genome::Model::Build->get($self->somatic_validation_build_id);
            $build_or_model ||= Genome::Model->get($self->somatic_validation_model_id); #fallback for now

            $build_or_model or die "Could not find model ($self->somatic_validation_model_id)\n";

            my $ref_seq_build_id = $build_or_model->reference_sequence_build->build_id;
            $ref_seq_build = Genome::Model::Build->get($ref_seq_build_id);
        } else {
            my $model = Genome::Model->get($self->tumor_val_model_id);
            my $ref_seq_build_id = $model->reference_sequence_build->build_id;
            $ref_seq_build = Genome::Model::Build->get($ref_seq_build_id);
        }
    }

    #create reference sequence using new contigs (define new reference and track new reference build)
    my $ref_model_id = $ref_seq_build->model_id;

    my $version = "500bp_assembled_contigs_sv";
    #don't overwrite an existing model...

    $version = checkRefBuildName($sample_id . "_SV_Contigs",$version);

    my $new_ref_cmd = Genome::Model::Command::Define::ImportedReferenceSequence->create(
        species_name => 'human',
        use_default_sequence_uri => '1',
        derived_from => $ref_seq_build,
        append_to => $ref_seq_build,
        version => $version,
        fasta_file => $contigs_file,
        prefix => $sample_id . "_SV_Contigs",
    );
    unless ($new_ref_cmd->execute) {
        $self->error_message('Failed to execute the definition of the new reference sequence with added contigs.');
        return;
    }
    my $new_ref_build_id = $new_ref_cmd->result_build_id;
    my $new_ref_build = Genome::Model::Build->get($new_ref_build_id);
    my $new_ref_event = $new_ref_build->the_master_event;
    my $new_ref_event_id = $new_ref_event->id;
    my $new_ref_event_class = $new_ref_event->class;
    while ($new_ref_event->event_status eq 'Running' || $new_ref_event->event_status eq 'Scheduled') {
        sleep 600;
        $new_ref_event = $new_ref_event_class->load($new_ref_event_id);
    }
    unless ($new_ref_event->event_status eq 'Succeeded') {
        $self->error_message('New reference build not successful.');
        return;
    }


    my $new_tumor_model_id;
    my $new_normal_model_id;


    if($input_model_type eq "somval"){
        #create new tumor and normal validation models to align data to new reference

        #first, grab a bunch of params we need:

        my ($build, $model);

        if($self->somatic_validation_build_id) {
            $build = Genome::Model::Build->get($self->somatic_validation_build_id);
            $model = $build->model;
        } else {
            $model = Genome::Model->get($self->somatic_validation_model_id);
            $build = $model->last_succeeded_build;
        }
        $model or
            die "Could not find model ($self->somatic_validation_model_id\n";

        $build or
            die "Could not find last succeeded build from somatic model $self->somatic_validation_model_id.\n";

        my $tumor_sample = $build->tumor_sample;
        my $tumor_sample_id = $build->tumor_sample->sample_id;
        my $tumor_subject = Genome::Subject->get($tumor_sample->subject_id);

        my $normal_sample = $build->normal_sample;
        my $normal_sample_id = $build->normal_sample->sample_id;
        my $normal_subject = Genome::Subject->get($normal_sample->subject_id);
        #instrument data ids, sample ids
        my @tumor_instrument_data;
        my @normal_instrument_data;
        my @instrument_data_ids = $build->instrument_data_ids;
        foreach my $instid (@instrument_data_ids){
            my $inst_data = Genome::InstrumentData->get($instid);
            my $library_id = $inst_data->library_id;
            my $library = Genome::Library->get($library_id);
            my $sample_id = $library->sample_id;

            if($sample_id eq $tumor_sample_id){
                push(@tumor_instrument_data,$inst_data);
            } elsif ($sample_id eq $normal_sample_id){
                push(@normal_instrument_data,$inst_data);

            } else {
                die "sample id $sample_id from instrument data does not match sample id from either tumor (" . $build->tumor_sample->sample_id . ") or normal (" . $build->normal_sample->sample_id . ")\n";
            }
        }

        #make sure these model names aren't taken. If they are, add a digit to the end
        my $new_tumor_model_name = $sample_id . "-Tumor-SV-Validation-ManRevContigs";
        my $new_normal_model_name = $sample_id . "-Normal-SV-Validation-ManRevContigs";

        $new_tumor_model_name = checkModelName($new_tumor_model_name);
        $new_normal_model_name = checkModelName($new_normal_model_name);

        print STDERR "creating models $new_normal_model_name and $new_tumor_model_name\n";

        my $new_pp = Genome::ProcessingProfile->get("2599983");
        print STDERR "---\n";
        print STDERR join("\n",($new_ref_build,join("|",@tumor_instrument_data),$tumor_subject,$new_tumor_model_name,$new_pp)) . "\n";
        print STDERR "---\n";



        #new tumor model
        my $tumor_copy = Genome::Model::Command::Define::ReferenceAlignment->create(
            reference_sequence_build => $new_ref_build,
            auto_build_alignments => 0,
            auto_assign_inst_data => 0,
            instrument_data => \@tumor_instrument_data,
            subject => $tumor_subject,
            model_name => $new_tumor_model_name,
            processing_profile => $new_pp
            );
        $tumor_copy->dump_status_messages(1);
        $tumor_copy->execute or die "tumor define failed";
        $new_tumor_model_id = $tumor_copy->result_model_id;
        my $new_tumor_model = Genome::Model->get($new_tumor_model_id);
        print STDERR "tumor model defined: $new_tumor_model_id\n";

        #new normal model
        my $normal_copy = Genome::Model::Command::Define::ReferenceAlignment->create(
            reference_sequence_build => $new_ref_build,
            auto_build_alignments => 0,
            auto_assign_inst_data => 0,
            instrument_data => \@normal_instrument_data,
            subject => $normal_subject,
            model_name => $new_normal_model_name,
            processing_profile => $new_pp
            );
        $normal_copy->dump_status_messages(1);
        $normal_copy->execute or die "copy failed";
        $new_normal_model_id = $normal_copy->result_model_id;
        my $new_normal_model = Genome::Model->get($new_normal_model_id);
        print STDERR "normal model defined: $new_normal_model_id\n";

    #-------------------------------------------------------
    } else { # $input_model_type eq "pairedref"

        #copy tumor and normal validation models to align data to new reference

        #make sure these model names aren't taken. If they are, add a digit to the end
        my $new_tumor_model_name = $sample_id . "-Tumor-SV-Validation-ManRevContigs";
        my $new_normal_model_name = $sample_id . "-Normal-SV-Validation-ManRevContigs";

        $new_tumor_model_name = checkModelName($new_tumor_model_name);
        $new_normal_model_name = checkModelName($new_normal_model_name);

        print STDERR "creating models $new_normal_model_name and $new_tumor_model_name\n";


        #my $new_pp = "dlarson bwa0.5.9 -q 5 indel contig test picard1.42";
        #my $new_pp = Genome::ProcessingProfile->get("2599983");
        my $new_pp = 2599983;
        my $tumor_model = Genome::Model->get($self->tumor_val_model_id) or die "Could not find tumor model with id $self->tumor_val_model_id.\n";
        my $normal_model = Genome::Model->get($self->normal_val_model_id) or die "Could not find normal model with id $self->normal_val_model_id.\n";

        #new tumor model
        my $tumor_copy = Genome::Model::Command::Copy->create(
            model => $tumor_model,
            overrides => [
                'name='.$new_tumor_model_name,
                'auto_build_alignments=0',
                'processing_profile='.$new_pp,
                'reference_sequence_build='.$new_ref_build_id,
                'annotation_reference_build=',
                'region_of_interest_set_name=',
                'dbsnp_build=',
            ],
            );
        $tumor_copy->dump_status_messages(1);
        $tumor_copy->execute or die "copy failed";
        my $new_tumor_model = $tumor_copy->_new_model;
        my $new_tumor_model_id = $new_tumor_model->id;

        #new normal model
        my $normal_copy = Genome::Model::Command::Copy->create(
            model => $normal_model,
            overrides => [
                'name='.$new_normal_model_name,
                'auto_build_alignments=0',
                'processing_profile='.$new_pp,
                'reference_sequence_build='.$new_ref_build_id,
                'annotation_reference_build=',
                'region_of_interest_set_name=',
                'dbsnp_build=',
            ],
            );
        $normal_copy->dump_status_messages(1);
        $normal_copy->execute or die "copy failed";
        my $new_normal_model = $normal_copy->_new_model;
        my $new_normal_model_id = $new_normal_model->id;
    }

    $self->realignment_model_ids([$new_tumor_model_id, $new_normal_model_id]);

    #final notices to user
    print "\n\n\n###################### IMPORTANT INFORMATION ######################\n\n";

    #alert user to run builds for these copied models
    print "To start alignments against the new reference sequence which contains the SV contigs, please run this command from genome stable:\n\ngenome model build start $new_tumor_model_id $new_normal_model_id\n\n";

    return 1;

}

1;


sub checkModelName{
    my $name = shift;

    my @models = Genome::Model->get('name LIKE' => $name . '%');

    my $max=-1;
    for my $model (@models){
        #if we have models with a suffix already, store the highest suffix
        if ($model->name =~ /$name-(\d+)/){
            if($1 > $max){
                $max = $1;
            }
            #else if we have a match at all for this model name
        } elsif ($model->name =~ /$name/){
            $max = 0;
        }
    }
    if($max > -1){
        $max++;
        $name = $name . "-$max";
    }
    return $name;
}

sub checkRefBuildName{
    my $sample_id = shift;
    my $version = shift;

    #FIXME this should be all about a specific build type...
    my @builds = Genome::Model::Build->get("model.name LIKE" => $sample_id . "-human%");

    my $max=-1;
    foreach my $build (@builds){
        my $v = $build->version;
        #if we have models with a suffix already, store the highest suffix
        if ($v=~/$version-(\d+)/){
            if($1 > $max){
                $max = $1;
            }
            #else if we have a match at all for this model name
        } elsif ($v=~/$version/){
            $max = 0;
        }
    }
    if($max > -1){
        $max++;
        $version = $version . "-$max";
    }
    return $version;
}

1;
