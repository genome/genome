package Genome::InstrumentData::AlignmentResult::Bfast;

use strict;
use warnings;
use File::Basename;

use Bio::SeqIO;
 
use Genome;

class Genome::InstrumentData::AlignmentResult::Bfast {
    is => 'Genome::InstrumentData::AlignmentResult',
    has_constant => [
        aligner_name => { value => 'Bfast', is_param=>1 },
    ],
    has_transient_optional => [
         _Bfast_sam_cmd => { is=>'Text' }
    ]
};

sub required_arch_os { 'x86_64' }

# fill me in here with what compute resources you need.
sub required_rusage { 
    my $self = shift;
    my %params = $self->decomposed_aligner_params;
    my $max = 0;
    foreach my $cmd (keys %params) {
        if ($params{$cmd} =~ /-n\s+([0-9]+)/) {
            if ($1 > $max) {
                $max = $1;
            }
        }
    }
    if ($max > 0) { 
        return "-R 'select[model!=Opteron250 && type==LINUX64 && tmp>90000 && mem>24000] span[hosts=1] rusage[tmp=90000, mem=24000]' -M 24000000 -n $max";
    } else {
        print "Could not determine number of cores to use from params! Defaulting to 4.";
        return "-R 'select[model!=Opteron250 && type==LINUX64 && tmp>90000 && mem>24000] span[hosts=1] rusage[tmp=90000, mem=24000]' -M 24000000 -n 4";
    }
}

sub _run_aligner {
    my $self = shift;
    my @input_pathnames = @_;
    
    my $tmp_dir = $self->temp_scratch_directory;
    my $fq_file = "$tmp_dir/merged.fq";
    my $tmp_matches_file = "$tmp_dir/bfast.matches.bmf";
    my $tmp_aligned_file = "$tmp_dir/bfast.aligned.baf";
    my $tmp_reported_file = "$tmp_dir/bfast.reported.sam";
    my $staging_sam_file = "$tmp_dir/all_sequences.sam";
    
    # get refseq info
    my $reference_build = $self->reference_build;
    my $ref_file = $reference_build->full_consensus_path('fa');
    #my $ref_basename = File::Basename::dirname($ref_file);
    #my $reference_fasta_path = sprintf("%s/%s", $reference_build->data_directory, $ref_basename);

    # original hard coded path
    #my $ref_file = "/gscmnt/sata820/info/medseq/alignment-test/bfast/reference/all_sequences.fa";
    my @ref_files = (
        $ref_file,
        $ref_file.".nt.1.1.bif",
        $ref_file.".nt.2.1.bif",
        $ref_file.".nt.3.1.bif",
        $ref_file.".nt.4.1.bif",
        $ref_file.".nt.5.1.bif",
        $ref_file.".nt.6.1.bif",
        $ref_file.".nt.7.1.bif",
        $ref_file.".nt.8.1.bif",
        $ref_file.".nt.9.1.bif",
        $ref_file.".nt.10.1.bif",
        $ref_file.".nt.brg"
    );

    foreach my $bif (@ref_files) {
        unless (-e $bif) {
            $self->error_message("Index $bif does not exist. Please create Bfast indexes in the correct location, or implement Genome::InstrumentData::Alignment->resolve_reference_build() in this module.");
            $self->die($self->error_message);
        }
    }
     
    my $bfast_path = Genome::Model::Tools::Bfast->path_for_bfast_version($self->aligner_version);

    my %aligner_params = $self->decomposed_aligner_params;

    #### STEP 1: If there are paired end reads, merge them into $fq_file. otherwise,
    ####    set $fq_file to the one given pathname from @input_pathnames
    
    print "Merging $input_pathnames[0] and $input_pathnames[1] into $fq_file\n";

    if (scalar(@input_pathnames) == 2) {
        my $fastq_in1 = $input_pathnames[0];
        my $fastq_in2 = $input_pathnames[1];
        my $fastq_out = $fq_file;
         
        my $in_1_obj = Bio::SeqIO->new(
            -file   => $fastq_in1,
            -format => 'fastq',
        );

        my $in_2_obj = Bio::SeqIO->new(
            -file   => $fastq_in2,
            -format => 'fastq',
        );

        my $out_obj = Bio::SeqIO->new(
            -file   => ">$fastq_out",
            -format => 'fastq',
        );

        my $read_count = 0;

        while (1) {
            $read_count++;
            my $fastq_obj1 = $in_1_obj->next_seq || last;
            my $fastq_obj2 = $in_2_obj->next_seq || last;
            my $id1 = $fastq_obj1->id;
            my $id2 = $fastq_obj2->id;
            chop($id1);
            chop($id1);
            chop($id2);
            chop($id2);
            if ($id1 ne $id2) {
                $self->error_message("reads at position $read_count $id1 and $id2 have different ids");
                $self->die($self->error_message);
            }

            $fastq_obj1->id($id1);
            $fastq_obj2->id($id2);

            $out_obj->write_fastq($fastq_obj1);
            $out_obj->write_fastq($fastq_obj2);
        }

        unless (-s $fq_file) {
            $self->error_message("Unable to merge $input_pathnames[0] and $input_pathnames[1] into $fq_file");
            $self->die($self->error_message);
        }
        # Bio::SeqIO handles should automatically be destroyed when we leave this scope
    } elsif (scalar(@input_pathnames) == 1) {
        $fq_file = $input_pathnames[0];
    } else {
        $self->error_message("number of input pathnames to Bfast was not 1 or 2 (this probably shouldn't happen)");
        $self->die($self->error_message);
    }

    #### STEP 2: Search the indexes
    
    {
        # note the slash after tmp_dir. bfast dies without a trailing slash
        my $cmdline = $bfast_path . sprintf(' match %s -T %s -f %s -r %s > %s',
            $aligner_params{match_params}, "$tmp_dir/", $ref_file, $fq_file, $tmp_matches_file);

        Genome::Sys->shellcmd(
            cmd             => $cmdline,
            input_files     => [ @ref_files, $fq_file ],
            output_files    => [ $tmp_matches_file ],
            skip_if_output_is_present => 0,
        );

        unless (-s $tmp_matches_file) {
            $self->error_message("Unable to search the indexes. Match file $tmp_matches_file is zero length, so something went wrong.");
            $self->die($self->error_message);
        }
    }
    
    #### STEP 3: Perform local alignment
    
    {
        my $cmdline = $bfast_path. sprintf(' localalign %s -f %s -m %s > %s',
            $aligner_params{localalign_params}, $ref_file, $tmp_matches_file, $tmp_aligned_file);

        Genome::Sys->shellcmd(
            cmd             => $cmdline,
            input_files     => [ @ref_files, $tmp_matches_file ],
            output_files    => [ $tmp_aligned_file ],
            skip_if_output_is_present => 0,
        );

        unless (-s $tmp_aligned_file) {
            $self->error_message("Unable to align. Aligned file $tmp_aligned_file is zero length, so something went wrong.");
            $self->die($self->error_message);
        }
    }
    
    #### STEP 4: Filter alignments
    
    {
        my $cmdline = $bfast_path . sprintf(' postprocess %s -f %s -i %s > %s',
            $aligner_params{postprocess_params}, $ref_file, $tmp_aligned_file, $tmp_reported_file);

        Genome::Sys->shellcmd(
            cmd             => $cmdline,
            input_files     => [ @ref_files, $tmp_aligned_file ],
            output_files    => [ $tmp_reported_file ],
            skip_if_output_is_present => 0,
        );

        unless (-s $tmp_reported_file) {
            $self->error_message("Unable to filter alignments. Sam file $tmp_reported_file is zero length, so something went wrong.");
            $self->die($self->error_message);
        }

        # put your output file here, append to this file!
            #my $output_file = $self->temp_staging_directory . "/all_sequences.sam"

        # cowboy debugging
        #print $cmdline;
        #<STDIN>;
        die "Failed to process sam command line, error_message is ".$self->error_message unless $self->_filter_sam_output($tmp_reported_file, $staging_sam_file);

        #<STDIN>;
    }

    # TODO (iferguso) should do something to handle the log files

    return 1;
}

sub fillmd_for_sam {
    # it appears that bfast does put in MD string, so...
    return 0;
}

sub _filter_sam_output {
    my ($self, $aligned_sam_file, $all_sequences_sam_file) = @_;
    # some aligners output unaligned reads to a separate file. bfast appears to put them in the same file

    my $aligned_fh = IO::File->new( $aligned_sam_file );
    if ( !$aligned_fh ) {
            $self->error_message("Error opening bfast output sam file for reading $!");
            return;
    }
    $self->debug_message("Opened $aligned_sam_file");

    my $all_seq_fh = IO::File->new(">>$all_sequences_sam_file");
    if ( !$all_seq_fh ) {
        $self->error_message("Error opening all seq sam file for writing $!");
        return;
    }
    $self->debug_message("Opened $all_sequences_sam_file");
    
    while (<$aligned_fh>) {
        #write out the aligned map, excluding the default header- all lines starting with @.
        $all_seq_fh->print($_) unless $_ =~ /^@/;
    }

    $aligned_fh->close;
    $all_seq_fh->close;
    return 1;
}

sub decomposed_aligner_params {
    my $self = shift;
    # TODO this may be redundant considering the conditionals below
    my $params = $self->aligner_params || "-n 4:-n 4:-n 4";
    
    my @spar = split /\:/, $params;
    # TODO this could potentially be a problem if we don't want to, say, force 4 cores when not otherwise specified
    if ($spar[0] !~ /-n/) { $spar[0] .= "-n 4"; }
    if ($spar[1] !~ /-n/) { $spar[1] .= "-n 4"; }
    if ($spar[2] !~ /-n/) { $spar[2] .= "-n 4"; }

    return ('match_params' => $spar[0], 'localalign_params' => $spar[1], 'postprocess_params' => $spar[2]);
}

sub aligner_params_for_sam_header {
    my $self = shift;
    
    my %params = $self->decomposed_aligner_params;
    
    return "bfast match $params{match_params}; bfast localalign $params{localalign_params}; bfast postprocess $params{postprocess_params}";

    # for bwa this looks like "bwa aln -t4; bwa samse 12345'
}
