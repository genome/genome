package Genome::InstrumentData::AlignmentResult::BwamemStream;

use strict;
use warnings;

use Data::Dumper;
use File::Basename;
use File::Copy qw(move);
use Path::Class;
use Genome;
use Getopt::Long;

class Genome::InstrumentData::AlignmentResult::BwamemStream {
    is => 'Genome::InstrumentData::AlignmentResult::BwamemBase',
    has_constant => [
        aligner_name => { value => 'bwamem-stream', is_param=>1 },
    ],
};


sub required_memory_gb {
    return 42;
}

sub fillmd_for_sam {
    return 0;
}

sub required_rusage {
    my $class = shift;
    my %p = @_;
    my $instrument_data = delete $p{instrument_data};
    my $aligner_params  = delete $p{aligner_params};

    my $tmp_mb = $class->tmp_megabytes_estimated($instrument_data);
    my $mem_mb = 1024 * $class->required_memory_gb;
    my $cpus = 4;

    my $params_hash = $class->_param_string_to_hash($aligner_params);
    if (exists $params_hash->{t}) {
        $cpus = $params_hash->{t};
    }

    my $mem_kb = $mem_mb*1024;
    my $tmp_gb = $tmp_mb/1024;

    my $queue = $ENV{GENOME_LSF_QUEUE_ALIGNMENT_DEFAULT};
    $queue = $ENV{GENOME_LSF_QUEUE_ALIGNMENT_PROD} if (Genome::Config->should_use_alignment_pd);
    #$queue = "alignment-pd";

    my $host_groups;
    $host_groups = qx(bqueues -l $queue | grep ^HOSTS:);
    $host_groups =~ s/\/\s+/\ /;
    $host_groups =~ s/^HOSTS:\s+//;

    my $select  = "select[ncpus >= $cpus && mem >= $mem_mb && gtmp >= $tmp_gb] span[hosts=1]";
    my $rusage  = "rusage[mem=$mem_mb, gtmp=$tmp_gb]";
    my $options = "-M $mem_kb -n $cpus -q $queue";

    my $required_usage = "-R '$select $rusage' $options";

    my $select_check = "select[ncpus >= $cpus && maxmem >= $mem_mb && maxgtmp >= $tmp_gb] span[hosts=1]";
    my $select_cmd = "bhosts -R '$select_check' $host_groups | grep ^blade";

    my @selected_blades = qx($select_cmd);

    if (@selected_blades) {
        return $required_usage;
    } else {
        die $class->error_message("Failed to find hosts that meet resource requirements ($required_usage). [Looked with `$select_cmd`]");
    }
}

sub tmp_megabytes_estimated {
    my $class = shift || die;
    my $instrument_data = shift;

    my $default_megabytes = 90000;


    if (not defined $instrument_data) {
        return $default_megabytes;
    } elsif ($instrument_data->bam_path) {
        my $bam_path = $instrument_data->bam_path;

        my $scale_factor = 6.0;

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

sub _faidx_indexed_fasta {
    my $self = shift;
    return $self->reference_build->full_consensus_path('fa');
}

sub _run_aligner {
    my $self = shift;
    my @input_paths = @_;

    $self->_verify_params_and_inputs(@input_paths);

    # Verify the aligner and get params.

    # get temp dir
    my $tmp_dir  = $self->temp_scratch_directory;
    my $log_path = $tmp_dir . '/aligner.log';
    my $out_sam  = $self->scratch_sam_file_path;

    my $param_hash = $self->decomposed_aligner_params;
    my $num_threads = delete $param_hash->{t} || 1;
    my $header_extra = $self->_sam_header_extra;
    if (exists $header_extra->{RG}) {
        my $rg_str = $header_extra->{RG};
        $rg_str =~ s/\t/\\t/g;
        $param_hash->{R} = $rg_str;
    }
    delete $param_hash->{R};

    my $param_string = $self->_param_hash_to_string($param_hash);

    # Run mem
    $self->debug_message("Running bwa mem.");

    my $bam_output_path = $self->temp_staging_directory . "/all_sequences.bam";

    my %params =  (
        sam_header_path => $self->scratch_sam_file_path,
        aligner_log_path => $log_path,
        num_threads => $num_threads,
        bwa_version => $self->aligner_version,
        input_fastqs => \@input_paths,
        output_file => $bam_output_path,
        indexed_fasta => $self->_faidx_indexed_fasta,
        aligner_index_fasta => $self->_aligner_index_fasta,
        aligner_params => $param_string,
        );

    $params{samtools_version} = $self->samtools_version if defined $self->samtools_version;

    my $cmd = Genome::Model::Tools::Bwa::RunMem->create(%params);


    die $self->error_message("bwa mem failed!") unless
        $cmd->execute;

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

    return 1;
}

1;
