package Genome::InstrumentData::AlignmentResult::Bwamem;

use strict;
use warnings;
use File::Copy qw/move/;
use Path::Class;
use Genome;
use Getopt::Long;

class Genome::InstrumentData::AlignmentResult::Bwamem {
    is => 'Genome::InstrumentData::AlignmentResult::BwamemBase',
    has_constant => [
        aligner_name => { value => 'bwamem', is_param=>1 },
    ],
};

sub required_memory_gb {
    return 16;
}

sub fillmd_for_sam {
    return 1;
}

sub required_rusage {
    my $class = shift;
    my %p = @_;
    my $instrument_data = delete $p{instrument_data};
    my $aligner_params  = delete $p{aligner_params};

    my $tmp_mb = $class->tmp_megabytes_estimated($instrument_data);
    my $mem_mb = 1024 * $class->required_memory_gb;
    my $cpus = 4;

    if ($aligner_params and $aligner_params =~ /-t\s*([0-9]+)/) {
        $cpus = $1;
    }

    my $mem_kb = $mem_mb*1024;
    my $tmp_gb = $tmp_mb/1024;

    my $queue = $ENV{GENOME_LSF_QUEUE_ALIGNMENT_DEFAULT};
    $queue = $ENV{GENOME_LSF_QUEUE_ALIGNMENT_PROD} if (Genome::Config->should_use_alignment_pd);

    my $host_groups;
    $host_groups = qx(bqueues -l $queue | grep ^HOSTS:);
    $host_groups =~ s/\/\s+/\ /;
    $host_groups =~ s/^HOSTS:\s+//;

    my $select  = "select[ncpus >= $cpus && mem >= $mem_mb && gtmp >= $tmp_gb] span[hosts=1]";
    my $rusage  = "rusage[mem=$mem_mb, gtmp=$tmp_gb]";
    my $options = "-M $mem_kb -n $cpus -q $queue";

    my $required_usage = "-R '$select $rusage' $options";

    #check to see if our resource requests are feasible (This uses "maxmem" to check theoretical availability)
    #factor of four is based on current six jobs per host policy this should be revisited later
    my $select_check = "select[ncpus >= $cpus && maxmem >= " . ($mem_mb * 4) . " && maxgtmp >= $tmp_gb] span[hosts=1]";
    my $select_cmd = "bhosts -R '$select_check' $host_groups | grep ^blade";

    my @selected_blades = qx($select_cmd);

    if (@selected_blades) {
        return $required_usage;
    } else {
        die $class->error_message("Failed to find hosts that meet resource requirements ($required_usage). [Looked with `$select_cmd`]");
    }
}

# TODO copied verbatim from normal bwa, but this may be totally off for bwa mem
sub tmp_megabytes_estimated {
    my $class = shift || die;
    my $instrument_data = shift;

    my $default_megabytes = 90000;


    if (not defined $instrument_data) {
        return $default_megabytes;
    } elsif ($instrument_data->bam_path) {
        my $bam_path = $instrument_data->bam_path;

        my $scale_factor = 3.25; # assumption: up to 3x during sort/fixmate/sort and also during fastq extraction (2x) + bam = 3
        # 3.25 scale factor worked for bwa but seems inefficient for bwamem;
        # we've seen usage up to 400 GB when 106 GB was requested.
        $scale_factor *= 4;

        my $bam_bytes = -s $bam_path;
        unless ($bam_bytes) {
            die $class->error_message("Instrument Data " . $instrument_data->id  . " has BAM ($bam_path) but has no size!");
        }

        if ($instrument_data->can('get_segments')) {
            my $bam_segments = scalar $instrument_data->get_segments;
            if ($bam_segments > 1) {
                $scale_factor = $scale_factor / $bam_segments;
            }
        }

        return int(($bam_bytes * $scale_factor) / 1024**2);
    } elsif ($instrument_data->can("calculate_alignment_estimated_kb_usage")) {
        my $kb_usage = $instrument_data->calculate_alignment_estimated_kb_usage;
        return int(($kb_usage * 3) / 1024) + 100; # assumption: 2x for the quality conversion, one gets rm'di after; aligner generates 0.5 (1.5 total now); rm orig; sort and merge maybe 2-2.5
    } else {
        return $default_megabytes;
    }

    return;
}

sub _run_aligner {
    my $self = shift;
    my @input_paths = @_;

    $self->_verify_params_and_inputs(@input_paths);

    # get temp dir
    my $tmp_dir  = $self->temp_scratch_directory;
    my $log_path = $tmp_dir . '/aligner.log';
    my $out_sam  = $self->scratch_sam_file_path;

    # Verify the aligner and get params.
    my $aligner_version = $self->aligner_version;
    my $cmd_path = Genome::Model::Tools::Bwa->path_for_bwa_version($aligner_version);
    my $param_hash = $self->decomposed_aligner_params;
    my $param_string = $self->_param_hash_to_string($param_hash);

    # Run mem
    $self->debug_message("Running bwa mem.");

    my $full_command = sprintf '%s mem %s %s %s 2>> %s',
        $cmd_path, $param_string, $self->_aligner_index_fasta,
        (join ' ', @input_paths), $log_path;
    $self->_stream_bwamem($full_command, $out_sam);

    # Verify the bwa mem logfile.
    unless ($self->_verify_bwa_mem_did_happen($log_path)) {
        die $self->error_message(
            "Error running bwa mem. Unable to verify a successful " .
            "run of bwa mem in the aligner log.");
    }

    # clean up the FASTQs in /tmp
    $self->debug_message("bwa mem command finished");
    $self->debug_message("Removing input FASTQs in tmp scratch space");
    $self->show_temporary_input_files_queue();
    $self->clear_temporary_input_files_queue();

    # Sort all_sequences.sam.
    $self->debug_message("Resorting all_sequences.sam by coordinate.");
    $self->_sort_sam($out_sam);

    return 1;
}

# Run bwa mem and stream through AddReadGroupTag
sub _stream_bwamem {
    my ($self, $full_command, $all_sequences) = @_;

    # Open pipe
    $self->debug_message("RUN: $full_command");
    $self->debug_message("Opening filehandle to stream output.");

    my $bwamem_fh = IO::File->new("$full_command |");

    unless ($bwamem_fh) {
        die $self->error_message(
            "Error running bwa mem. Unable to open filehandle " .
            "to stream bwa mem output.");
    }

    # Add RG tags
    $self->debug_message("Starting AddReadGroupTag.");
    my $all_sequences_fh = IO::File->new(">> $all_sequences");

    unless ($all_sequences_fh) {
        die $self->error_message(
            "Error running bwa mem. Unable to open all_sequences.sam " .
            "filehandle for AddReadGroupTag.");
    }

    # Prepare the RG command
    my $add_rg_cmd = Genome::Model::Tools::Sam::AddReadGroupTag->create(
       input_filehandle  => $bwamem_fh,
       output_filehandle => $all_sequences_fh,
       read_group_tag    => $self->read_and_platform_group_tag_id,
       pass_sam_headers  => 0,
    );

    # Disconnect from db
    $self->_disconnect_from_db();

    # Run RG command
    unless ($add_rg_cmd->execute) {
        die $self->error_message(
            "Error running bwa mem. AddReadGroupTag failed to execute.");
    }

    # Hopefully we're still disconnected
    $self->_check_db_connection();

    # Clean up
    $all_sequences_fh->close();
    $bwamem_fh->close();
    my $rv = $?;
    my $exit_code = $rv >> 8;
    my $core_dump = $rv & 128;

    if ($exit_code) {
        die $self->error_message(
            "Error running bwa mem. Expected exit code of '0' " .
            "but got $exit_code instead (\$? set to $rv).");
    }
    if ($core_dump) {
        die $self->error_message(
            "Error running bwa mem. Detected a coredump from " .
            "bwa mem (\$? set to $rv).");
    }
}

# Sort a sam file.
sub _sort_sam {
    my ($self, $given_sam) = @_;

    my $unsorted_sam = "$given_sam.unsorted";

    # Prepare sort command
    unless (move($given_sam, $unsorted_sam)) {
        die $self->error_message(
            "Unable to move $given_sam to $unsorted_sam. " .
            "Cannot proceed with sorting.");
    }

    # Memory is set to 80% of our LSF request; GMT::Picard should handle
    # non-integer memory requests correctly
    my $picard_sort_cmd = Genome::Model::Tools::Picard::SortSam->create(
        sort_order             => 'coordinate',
        input_file             => $unsorted_sam,
        output_file            => $given_sam,
        max_records_in_ram     => 2000000,
        maximum_memory         => ($self->required_memory_gb*0.8),
        maximum_permgen_memory => 256,
        temp_directory         => $self->temp_scratch_directory,
        use_version            => $self->picard_version,
    );

    # Disconnect from db
    $self->_disconnect_from_db();

    # Run sort command
    unless ($picard_sort_cmd and $picard_sort_cmd->execute) {
        die $self->error_message(
            "Failed to create or execute Picard sort command.");
    }

    # Hopefully we're still disconnected
    $self->_check_db_connection();

    # Clean up
    unless (unlink($unsorted_sam)) {
        $self->debug_message("Could not unlink $unsorted_sam.");
    }

    return $given_sam;
}

1;
