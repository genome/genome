package Genome::InstrumentData::AlignmentResult::PerLaneTophat;

use strict;
use warnings;

use Genome;

class Genome::InstrumentData::AlignmentResult::PerLaneTophat {
    is => 'Genome::InstrumentData::AlignmentResult',
    has_constant => [
        aligner_name => {
            value => 'tophat',
            is_param => 1
        },
    ],
    has_calculated => [
        bowtie_version => {
            is => "Text",
            calculate => \&_get_bowtie_version,
            calculate_from => ['class', 'aligner_params'],
            doc => 'the version of bowtie passed in - will need to be translated into appropriate tophat params'
        },
    ],
};

sub required_arch_os { 'x86_64' }

#WARNING: THIS RESOURCE REQUEST IS ONLY VALID FOR EVENT-BASED WORKFLOWS.
#IF YOU WANT TO CHANGE THE RESOURCE REQUEST FOR NORMAL WORKFLOWS, EDIT
#THE FILE Genome::InstrumentData::Command::AlignReads::PerLaneTophat INSTEAD!!!
sub required_rusage {
    my $class = shift;
    
    my $mem_mb = 32000;
    my $mem_kb = $mem_mb*1024;
    my $tmp_gb = 400;
    my $cpus = 4;
    
    my $select  = "select[ncpus >= $cpus && mem >= $mem_mb && gtmp >= $tmp_gb] span[hosts=1]";
    my $rusage  = "rusage[mem=$mem_mb, gtmp=$tmp_gb]";
    my $options = "-M $mem_kb -n $cpus";
    my $required_usage = "-R '$select $rusage' $options";
    
    return $required_usage;
}

sub _run_aligner {
    my $self = shift;
    my @input_pathnames = @_;


    # get refseq info
    my $reference_build = $self->reference_build;
    # This is your scratch directory.  Whatever you put here will be wiped when the alignment
    # job exits.
    my $scratch_directory = $self->temp_scratch_directory;
    # This is the alignment output directory.  Whatever you put here will be synced up to the
    # final alignment directory that gets a disk allocation.
    my $staging_directory = $self->temp_staging_directory;

    my $tophat_cmd = $self->_get_tophat_cmd(\@input_pathnames);

    # disconnect the db handle before this long-running event
    if (Genome::DataSource::GMSchema->has_default_handle) {
        $self->status_message("Disconnecting GMSchema default handle.");
        Genome::DataSource::GMSchema->disconnect_default_dbh();
    }

    Genome::Sys->shellcmd(
        cmd => $tophat_cmd,
        input_files => \@input_pathnames,
        output_files => [ "$staging_directory/accepted_hits.bam", "$staging_directory/unmapped.bam" ]
    );

    rename("$staging_directory/accepted_hits.bam", "$scratch_directory/accepted_hits.bam");
    rename("$staging_directory/unmapped.bam", "$scratch_directory/unmapped.bam");

    my $bam_with_unaligned_reads_cmd = Genome::Model::Tools::Picard::MergeSamFiles->create(
        input_files => ["$scratch_directory/unmapped.bam", "$scratch_directory/accepted_hits.bam",],
        output_file => "$scratch_directory/all_reads.bam",
        use_version => $self->picard_version,
    );

    unless($bam_with_unaligned_reads_cmd->execute()){
        die($self->error_message("Unable to create a merged bam with both aligned and unaligned reads."));
    }

    my $sam_path = Genome::Model::Tools::Sam->path_for_samtools_version($self->samtools_version);
    my $sam_cmd = "$sam_path view $scratch_directory/all_reads.bam |";
    my $add_rg_cmd = Genome::Model::Tools::Sam::AddReadGroupTag->create(
      input_filehandle => IO::File->new($sam_cmd),
      output_filehandle => $self->_sam_output_fh,
      read_group_tag => $self->read_and_platform_group_tag_id,
      pass_sam_headers => 0,
    );

    unless($add_rg_cmd->execute){
      $self->error_message("Adding read group to sam file failed!");
      die $self->error_message;
    }

    #promote other misc tophat result files - converted sam will be handled downstream
    rename("$scratch_directory/junctions.bed", "$staging_directory/junctions.bed");
    rename("$scratch_directory/insertions.bed", "$staging_directory/insertions.bed");
    rename("$scratch_directory/deletions.bed", "$staging_directory/deletions.bed");
    rename("$scratch_directory/logs", "$staging_directory/logs");
    return 1;
}

sub _check_read_count {
    my ($self) = @_;
    my $fq_rd_ct = $self->_fastq_read_count;
    my $sam_path = Genome::Model::Tools::Sam->path_for_samtools_version($self->samtools_version);

    my $cmd = "$sam_path view -F 256 -c " . $self->temp_staging_directory . "/all_sequences.bam";
    my $bam_read_count = `$cmd`;
    my $check = "Read count from bam: $bam_read_count and fastq: $fq_rd_ct";

    unless ($fq_rd_ct == $bam_read_count) {
        $self->error_message("$check does not match.");
        return;
    }
    $self->status_message("$check matches.");
    return 1;
}

sub alignment_bam_file_paths {
    my $self = shift;

    return ($self->output_dir . "/all_sequences.bam");
}

sub aligner_params_for_sam_header {
    my $self = shift;
    return "tophat " . $self->aligner_params;
}

sub aligner_params_required_for_index {
    return 1;
}

# Does your aligner set MD tags?  If so this should return 0, otherwise 1
sub fillmd_for_sam {
    return 0;
}

# If your aligner adds read group tags, or you are handling it in the wrapper, this needs to be 0
# otherwise 1.  If you are streaming output straight to BAM then you need to take care of adding RG
# tags with either the wrapper or the aligner itself, and this needs to be 0.
sub requires_read_group_addition {
    return 0;
}

# if you are streaming to bam, set this to 1.  Beware of read groups.
sub supports_streaming_to_bam {
    return 1;
}

# If your aligner accepts BAM files as inputs, return 1.  You'll get a set of BAM files as input, with
# suffixes to define whether it's paired or unparied.
# [input_file.bam:0] implies SE, [input_file.bam:1, input_file.bam:2] implies paired.
sub accepts_bam_input {
    return 0;
}

sub prepare_annotation_index {
    my $class = shift;
    my $annotation_index = shift;
    
    my $refindex = Genome::Model::Build::ReferenceSequence::AlignerIndex->get_with_lock(
        aligner_name => $annotation_index->aligner_name,
        aligner_version => $annotation_index->aligner_version,
        aligner_params => $annotation_index->aligner_params,
        reference_build => $annotation_index->reference_build,
    );

    unless ($refindex) {
        $class->error_message('Failed to find reference index!');
        return;
    }
    
    my $bowtie_extension = 'bt2';
    my $reference_extension = 'fa';
    my $aligner_params = $annotation_index->aligner_params;
    my $bowtie_version = $class->_get_bowtie_version($aligner_params);
    $aligner_params =~ s/--bowtie-version(?:\s+|=)(.+?)(?:\s+|$)//i;

    unless ($bowtie_version =~ /^2/) {
        $aligner_params .= ' --bowtie1';
        $bowtie_extension = 'ebwt';
        $reference_extension = 'bowtie';
    }

    my $reference_prefix = $refindex->full_consensus_path($reference_extension);
    my $staging_dir = $annotation_index->temp_staging_directory;
    
    # Build a transcriptome index if an annotation build is provided
    my $annotation_build = $annotation_index->annotation_build;

    my $gtf_path = $annotation_build->annotation_file('gtf',$refindex->reference_build_id);

    unless (defined($gtf_path)) {
        $class->error_message('There is no annotation GTF file defined for annotation_reference_transcripts build: '. $annotation_build->__display_name__);
        return;
    }
    
    unless (-s $gtf_path) {
        $class->error_message(sprintf("Annotation file %s does not exist", $gtf_path));
        return;
    }

    $class->status_message(sprintf("Confirmed non-zero annotation file is %s", $gtf_path));
    unless (symlink($gtf_path, sprintf("%s/all_sequences.gtf", $staging_dir))) {
        $class->error_message("Couldn't symlink annotation into the staging directory");
        return;
    }

    if ($aligner_params =~ /-G/) {
        $class->error_message('Annotation build \''. $annotation_build->__display_name__ .'\' is defined, but there seems to be a GTF file already in the read_aligner_params: '. $aligner_params);
        return;
    }
    if (defined($aligner_params)) {
        $aligner_params .= ' -G '. $gtf_path;
    } else {
        $aligner_params = ' -G '. $gtf_path;
    }

    my $transcriptome_index_prefix = $staging_dir .'/all_sequences';
    $aligner_params .= ' --transcriptome-index '. $transcriptome_index_prefix;

    my $path = Genome::Model::Tools::Tophat->path_for_tophat_version($refindex->aligner_version);
    my $bowtie_path = Genome::Model::Tools::Bowtie->path_variable_for_bowtie_version($bowtie_version);

    my $temp_directory = Genome::Sys->create_temp_directory();

    #a few in silico reads for gene TTN
    my $reference_fasta = $refindex->full_consensus_path('fa');
    $class->status_message("Generating fake reads from $reference_fasta and $gtf_path");
    my @fake_reads = $class->generate_fake_reads($reference_fasta, $gtf_path, 100, 50);
    my $fake_read_index = 1;
    @fake_reads = map { sprintf ">read%d\n%s\n", $fake_read_index++, $_ } @fake_reads;

    my ($fake_read) = Genome::Sys->create_temp_file_path('read.fasta');
    Genome::Sys->write_file($fake_read, @fake_reads);

    my $cmd = 'PATH=' . $bowtie_path . ':$PATH ' . $path .' '. $aligner_params .' --transcriptome-only --output-dir '. $temp_directory .' '. $reference_prefix .' '. $fake_read;

    my @output_files = (
            $transcriptome_index_prefix .'.gff',
            $transcriptome_index_prefix .'.fa',
            $transcriptome_index_prefix .'.1.'. $bowtie_extension,
            $transcriptome_index_prefix .'.2.'. $bowtie_extension,
            $transcriptome_index_prefix .'.3.'. $bowtie_extension,
            $transcriptome_index_prefix .'.4.'. $bowtie_extension,
            $transcriptome_index_prefix .'.rev.1.'. $bowtie_extension,
            $transcriptome_index_prefix .'.rev.2.'. $bowtie_extension,
    );

    Genome::Sys->shellcmd(
        cmd => $cmd,
        input_files => [$fake_read,$gtf_path],
        output_files => \@output_files,
    );
    return 1;
}


# This is used to generate fake reads given a reference fasta and a transcript gtf
sub generate_fake_reads {
    my ($class, $reference_fasta, $gtf_path, $read_length, $read_count) = @_;

    unless (-s $reference_fasta) {
        $class->error_message("The file $reference_fasta is empty or does not exist.");
        return;
    }

    unless (-s $gtf_path) {
        $class->error_message("The file $gtf_path is empty or does not exist.");
        return;
    }

    my $fa_fh = new Bio::DB::Fasta($reference_fasta);
    my @ids = $fa_fh->ids();

    open my $gtf_fh, '<', $gtf_path;

    my @reads;

    my $buffer = 5;
    my $seqid_limit = 5;

    # keep track of which chromosomes we've made reads for
    my %seqids_used;

    while (my $line = <$gtf_fh>) {
        last if @reads >= $read_count;

        chomp $line;

        my @fields = split /\t/, $line;

        my $seqid = $fields[0];
        my $kind  = $fields[2];
        my $start = $fields[3];
        my $stop  = $fields[4];

        next if (exists $seqids_used{$seqid}) and ($seqids_used{$seqid} >= $seqid_limit);
        next unless grep { $seqid eq $_ } @ids;
        next unless $kind eq 'CDS' or $kind eq 'exon';
        next if $stop - $start < $read_length + (3 * $buffer);

        my $seq = $fa_fh->seq($seqid, $start + $buffer, $start + $read_length + $buffer - 1);

        next if $seq =~ /[^ACTG]/;
        next if grep { $seq eq $_ } @reads;

        $seqids_used{$seqid}++;

        push @reads, $seq;
    }

    close $gtf_fh;

    # hard coded reads we were using before, in case our generated fake reads are no good
    push @reads, (
        "GATGTCTTTCAGCATGGAAGTAATCATGATTGGGGTGACTGTCAGGGTGGCACTGACTTGGTCATTGCCACAGACAAATGTGTATTCTGCCGAGTCCTCT",
        "GCAGCTGGTGTGTCAAGCACTTTAACACTCACAGTTGAAGACTTAGGTTCACCAACTCCATTTTCTATTGTGAGAATATATTTTCCTGTATCATTGCGAG",
        "TCCTTCTGTGCATGAGTGTTCTGAAGGGACTAGGGGCTCATAGTTTACCTGAGAGATCATGACATCAGGACTCTGGAGACTCTCCACGTGTCCCTCAGCT",
    );

    return @reads;
}

# Use this to prep the index.  Indexes are saved for each combo of aligner params & version, as runtime
# params for some aligners also call for specific params when building the index.
sub prepare_reference_sequence_index {
    my $class = shift;
    my $refindex = shift;

    # If you need the parameters the aligner is being run with, in order to customize the index, here they are.
    my $aligner_params = $refindex->aligner_params;

    my $staging_dir = $refindex->temp_staging_directory;
    my $staged_fasta_file = sprintf("%s/all_sequences.fa", $staging_dir);

    Genome::Sys->create_symlink($refindex->reference_build->get_sequence_dictionary("sam"), $staging_dir ."/all_sequences.dict" );
    my $bowtie_version = $class->_get_bowtie_version($aligner_params);
    my $bowtie_index = Genome::Model::Build::ReferenceSequence::AlignerIndex->get_or_create(
        reference_build_id => $refindex->reference_build_id,
        aligner_name => 'bowtie',
        aligner_version => $bowtie_version,
        test_name => $ENV{GENOME_ALIGNER_INDEX_TEST_NAME},
    );

    for my $filepath (glob($bowtie_index->output_dir . "/*")){
        my $filename = File::Basename::fileparse($filepath);
        next if $filename eq 'all_sequences.fa';
        next if $filename eq 'all_sequences.dict';
        Genome::Sys->create_symlink($filepath, $staging_dir . "/$filename");
    }

    $bowtie_index->add_user(
        label => 'uses',
        user => $refindex
    );

    my $actual_fasta_file = $staged_fasta_file;

    if (-l $staged_fasta_file) {
        $class->status_message(sprintf("Following symlink for fasta file %s", $staged_fasta_file));
        $actual_fasta_file = readlink($staged_fasta_file);
        unless($actual_fasta_file) {
            $class->error_message("Can't read target of symlink $staged_fasta_file");
            return;
        }

        if ($bowtie_version =~ /^2/) {
            my $bowtie_fasta_path = $bowtie_index->full_consensus_path('fa');
            Genome::Sys->create_symlink($bowtie_fasta_path, $staged_fasta_file . ".fa");
        } else {
            my $bowtie_fasta_path = $bowtie_index->full_consensus_path('bowtie');
            unless(-e $bowtie_fasta_path){
              $bowtie_fasta_path = $bowtie_index->full_consensus_path('fa');
            }
            Genome::Sys->create_symlink($bowtie_fasta_path, $staging_dir .'/all_sequences.bowtie.fa');
        }
    }

    return 1;
}

sub _get_bowtie_version {
    my ($class, $aligner_params) = @_;
    $aligner_params =~ /--bowtie-version(?:\s+|=)(.+?)(?:\s+|$)/i;
    return $1;
}

sub get_annotation_index {
    my $self = shift;
    my $annotation_build = shift || $self->annotation_build;

    my $index = Genome::Model::Build::ReferenceSequence::AnnotationIndex->get_with_lock(
        aligner_name => $self->aligner_name,
        aligner_version => $self->aligner_version,
        aligner_params => $self->aligner_params,
        reference_build => $self->reference_build,
        annotation_build => $annotation_build,
    );

    if (!$index) {
        die $self->error_message(sprintf("No annotation index prepared for %s with params %s and reference build %s with annotation build %s", $self->aligner_name, $self->aligner_params, $self->reference_build->id,$self->annotation_build->id));
    }

    return $index;
}

sub _get_modified_tophat_params {
    my $self = shift;
    my $params = $self->aligner_params;
    $params =~ s/--bowtie-version(?:\s+|=)(.+?)(?:\s+|$)//i;
    unless ($1 =~ /^2/) {
        $params .= " --bowtie1";
    }

    if ($self->annotation_build) {
        #TODOthis method doesn't exist
        my $annotation_index = $self->get_annotation_index;
        if ($params =~ /-G/) {
            die ('This processing_profile is requesting annotation_reference_transcripts \''. $self->annotation_build->__display_name__ .'\', but there seems to be a GTF file already defined in the read_aligner_params: '. $params);
        }
        #TODO: get ref index with annotation build
        my $transcriptome_index_prefix = $annotation_index->full_consensus_path();
        $params .= ' --transcriptome-index '. $transcriptome_index_prefix;
    }

    my $instrument_data = $self->instrument_data;
    # This is hard-coded for paired-end reads
    # TODO: Add support for fragment alignments
    # Default of 300 should be reasonable for 500bp PE libraries sequenced 2x100bp
    my $estimated_library_size = 350;
    my $estimated_sd = 50;
    
    #my $median_inner_insert_size = 300;
    #if ($instrument_data->resolve_median_insert_size) {
    #    $median_inner_insert_size  = ($instrument_data->resolve_median_insert_size - ($instrument_data->read_length * 2) );
    #}
    
    my $median_inner_insert_size = ($estimated_library_size - ($instrument_data->read_length * 2) );
    
    #my $sd_insert_size = 20;
    my $sd_insert_size = $estimated_sd;
    #if ($instrument_data->resolve_sd_insert_size) {
    #    $sd_insert_size = $instrument_data->resolve_sd_insert_size;
    #}
    $params .= ' --mate-inner-dist '. $median_inner_insert_size;
    $params .= ' --mate-std-dev '. $sd_insert_size;
    return $params;
}

sub _get_tophat_cmd {
    my $self = shift;
    my $input_filenames = shift;
    die("More than 2 input files cannot be specified!") unless $#$input_filenames < 2;

    my $path = Genome::Model::Tools::Tophat->path_for_tophat_version($self->aligner_version);
    my $bowtie_path = Genome::Model::Tools::Bowtie->path_variable_for_bowtie_version($self->bowtie_version); #tophat uses `which` to find executable--place desired version first
    my $params = $self->_get_modified_tophat_params;
    my $cmd = 'PATH=' . $bowtie_path . ':$PATH ' . $path . " " . $params . " -o " . $self->temp_staging_directory;

    my $bowtie_prefix = 'fa';
    if ($params =~ /--bowtie1/) {
        $bowtie_prefix = 'bowtie';
    }

    $cmd .=  " "  . $self->get_reference_sequence_index->full_consensus_path($bowtie_prefix) . " " . join(" " , @$input_filenames);
    return $cmd;
}

sub _use_alignment_summary_cpp { return; };

1;
