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
    my ($class, $annotation_index) = @_;
    unless (defined($annotation_index)) {
        $class->error_message("Required argument annotation_index was not passed!");
        return;
    }

    my $gtf_file = $class->_get_gtf_file($annotation_index) or return;
    my $reference_fasta = $class->_get_reference_fasta($annotation_index) or return;
    my $fake_reads_file = $class->_generate_fake_reads_file($reference_fasta, $gtf_file) or return;
    $class->_run_tophat_to_generate_annotation_index($annotation_index, $fake_reads_file, $reference_fasta, $gtf_file) or return;

    return 1;
}

sub _get_reference_fasta {
    my ($class, $annotation_index) = @_;

    # get the reference_fasta
    my $reference_index = Genome::Model::Build::ReferenceSequence::AlignerIndex->get_with_lock(
        aligner_name => $annotation_index->aligner_name,
        aligner_version => $annotation_index->aligner_version,
        aligner_params => $annotation_index->aligner_params,
        reference_build => $annotation_index->reference_build,
    );
    unless ($reference_index) {
        $class->error_message('Failed to find the reference index to retrieve the reference FASTA file!  '.
            "aligner_name: ".$annotation_index->aligner_name." aligner_version: ".$annotation_index->aligner_version.
            "aligner_params: ".$annotation_index->aligner_params." reference_build: ".$annotation_index->reference_build->id);
        return;
    }
    my $reference_fasta = $reference_index->full_consensus_path('fa');

    # complain if we can't get it
    my $msg = sprintf(
        "reference fasta (%s) from reference index (%s)",
        $reference_fasta, $reference_index->id,
    );
    if (-s $reference_fasta) {
        $class->debug_message("Found " .$msg);
    } else {
        $class->error_message("Failed to find " . $msg);
        return;
    }
    return $reference_fasta;
}

sub _generate_fake_reads_file {
    my ($class, $reference_fasta, $gtf_file) = @_;

    # call tool to make fake reads
    my $tmp_dir = Genome::Sys->create_temp_directory();
    my $fake_reads_file = File::Spec->join($tmp_dir, 'fake_reads.fastq');
    my $cmd = Genome::Model::Tools::BioSamtools::SimulateRnaSeqReads->create(
        reference_fasta_file => $reference_fasta,
        annotation_gtf_file => $gtf_file,
        fastq_file => $fake_reads_file,
        paired_end => 0,
        approximate_depth => 3,
        num_transcripts => 3,
    );

    unless ($cmd->execute()) {
        $class->error_message("Failed to generate fake reads.");
        return;
    }

    unless (-s $fake_reads_file) {
        $class->error_message("Fake reads file (%s) shouldn't be empty.");
        return;
    }
    return $fake_reads_file;
}

sub _get_gtf_file {
    my ($class, $annotation_index) = @_;

    my $annotation_build = $annotation_index->annotation_build;
    my $gtf_file = $annotation_build->annotation_file('gtf', $annotation_index->reference_build->id);
    my $fallback_gtf_file = $annotation_build->annotation_file('gtf');

    # complain if we can't get it
    my $msg = sprintf(
        "gtf file (%s) from annotation_build (%s) using reference id (%s)",
        $gtf_file, $annotation_build->id, $annotation_index->reference_build->id,
    );
    if (-s $gtf_file) {
        $class->debug_message("Found " .$msg);
        return $gtf_file;
    } elsif (-s $fallback_gtf_file) {
        $class->warning_message('Did not find %s. Falling back to %s', $msg, $fallback_gtf_file);
        return $fallback_gtf_file;
    } else {
        $class->error_message("Failed to find " . $msg);
        return;
    }
}

sub _run_tophat_to_generate_annotation_index {
    my ($class, $annotation_index, $fake_reads_file, $reference_fasta, $gtf_file) = @_;

    my $index_prefix = File::Spec->join($annotation_index->temp_staging_directory, 'all_sequences');
    my $aligner_params = $class->_get_aligner_params_to_generate_annotation_index($annotation_index, $gtf_file, $index_prefix) or return;

    # symlink in the gtf to live beside the index.
    Genome::Sys->create_symlink($gtf_file, $index_prefix . '.gtf');

    my $bowtie_path = $class->_get_bowtie_path_to_generate_annotation_index($annotation_index) or return;
    my $tophat_executable = Genome::Model::Tools::Tophat->path_for_tophat_version($annotation_index->aligner_version);

    my $cmd = "PATH=$bowtie_path" . ':$PATH' . " $tophat_executable $aligner_params $reference_fasta $fake_reads_file";

    my $bowtie_extension = 'bt2';
    my @output_files = (
            $index_prefix .'.gff',
            $index_prefix .'.fa',
            $index_prefix .'.1.'. $bowtie_extension,
            $index_prefix .'.2.'. $bowtie_extension,
            $index_prefix .'.3.'. $bowtie_extension,
            $index_prefix .'.4.'. $bowtie_extension,
            $index_prefix .'.rev.1.'. $bowtie_extension,
            $index_prefix .'.rev.2.'. $bowtie_extension,
    );

    Genome::Sys->shellcmd(
        cmd => $cmd,
        input_files => [$tophat_executable, $reference_fasta, $fake_reads_file],
        output_files => \@output_files,
    );
    return 1;
}

sub _get_bowtie_path_to_generate_annotation_index {
    my ($class, $annotation_index) = @_;

    my $aligner_params = $annotation_index->aligner_params;
    unless ($aligner_params) {
        $class->error_message("Aligner params cannot be empty (you must at the very least specify --bowtie-version)");
        return;
    }

    # determine the bowtie-version and remove it from aligner_params
    my $bowtie_version = $class->_get_bowtie_version($aligner_params);

    unless (defined($bowtie_version)) {
        $class->error_message("Couldn't find bowtie-version in aligner_params ($aligner_params)");
        return;
    }
    my $bowtie_path = Genome::Model::Tools::Bowtie->path_variable_for_bowtie_version($bowtie_version) or return;
    return $bowtie_path;
}

sub _aliger_params_has_gtf {
    my $params = shift;
    return ($params =~ /(^| )(-G|--GTF)/ ? 1 : 0);
}

sub _get_aligner_params_to_generate_annotation_index {
    my ($class, $annotation_index, $gtf_file, $index_prefix) = @_;

    my $aligner_params = $annotation_index->aligner_params;

    # remove --bowtie-version from aligner_params
    $aligner_params = $class->_remove_bowtie_version($aligner_params);
 
    # add in the options that tell TopHat that we want to store the index
    $aligner_params .= "--transcriptome-only";
    $aligner_params .= " --transcriptome-index '$index_prefix'";

    # add -G <gtf_file> option to tell TopHat what to make the index out of
    if (_aliger_params_has_gtf($aligner_params)) {
        my $annotation_build = $annotation_index->annotation_build;
        $class->error_message('Annotation build \''. $annotation_build->__display_name__ .
            '\' is defined, but there seems to be a GTF file already in the read_aligner_params: ' .
            $aligner_params);
        return;
    } else {
        $aligner_params .= ' -G '. $gtf_file;
    }

    # we con't care about the actual results of the 'fake' alignment so put them in a temp directory.
    my $temp_directory = Genome::Sys->create_temp_directory();
    $aligner_params .= " --output-dir '$temp_directory'";

    return $aligner_params;
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

our $BOWTIE_VERSION_REGEX = qr{--bowtie-version(?:\s+|=)(.+?)(?:\s+|$)};
sub _get_bowtie_version {
    my ($class, $aligner_params) = @_;
    $aligner_params =~ /$BOWTIE_VERSION_REGEX/i;
    return $1;
}

sub _remove_bowtie_version {
    my ($class, $aligner_params) = @_;
    $aligner_params =~ s/$BOWTIE_VERSION_REGEX//i;
    return $aligner_params;
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
    my $bowtie_version = $self->_get_bowtie_version($params);
    $params = $self->_remove_bowtie_version($params);
    unless ($bowtie_version =~ /^2/) {
        $params .= " --bowtie1";
    }

    if ($self->annotation_build) {
        #TODOthis method doesn't exist
        my $annotation_index = $self->get_annotation_index;
        if (_aliger_params_has_gtf($params)) {
            die ('This processing_profile is requesting annotation_reference_transcripts \''. $self->annotation_build->__display_name__ .'\', but there seems to be a GTF file already defined in the read_aligner_params: '. $params);
        }
        #TODO: get ref index with annotation build
        my $transcriptome_index_prefix = $annotation_index->full_consensus_path();
        $params .= ' --transcriptome-index '. $transcriptome_index_prefix;
    }

    #set number of threads automatically
    my $cpu_count = $self->_available_cpu_count;
    $self->status_message("CPU count is $cpu_count");
    if($params =~ /(^| )--num-threads[ =]\d+/) {
        $params =~ s/(^| )--num-threads[ =]\d+/$1--num-threads $cpu_count/;
    } elsif ($params =~ /(^| )-p[ =]\d+/) {
        $params =~ s/(^| )-p[ =]\d+/$1--num-threads $cpu_count/;
    } else {
        $params .= " --num-threads $cpu_count";
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
    if ($params =~ /(^| )--bowtie1/) {
        $bowtie_prefix = 'bowtie';
    }

    $cmd .=  " "  . $self->get_reference_sequence_index->full_consensus_path($bowtie_prefix) . " " . join(" " , @$input_filenames);
    return $cmd;
}

sub _use_alignment_summary_cpp { return; };

1;
