package Genome::InstrumentData::AlignmentResult::Merged::BamQc;

use strict;
use warnings;

use Genome;
use Sys::Hostname;
use File::Path;

class Genome::InstrumentData::AlignmentResult::Merged::BamQc {
    is => ['Genome::SoftwareResult::Stageable'],
    has_input => [
        alignment_result_id => {
            is  => 'Number',
            doc => 'ID of the result for the alignment data upon which to run coverage stats',
        },
    ],
    has_param => [
        picard_version => {
            is  => 'Text',
            doc => 'The version of Picard to use.',
        },
        samstat_version => {
            is  => 'Text',
            doc => 'The version of SamStat to use.',
        },
        fastqc_version => {
            is  => 'Text',
            doc => 'The version of FastQC to use.',
        },
        samtools_version => {
            is  => 'Text',
            doc => 'The version of Samtools to use.',
        },
        error_rate => {
            is  => 'Boolean',
            doc => 'Whether or not to generate error rate summaries.',
        },
        error_rate_version => {
            is  => 'Text',
            doc => 'The version of bam error rate C tool to use.',
        },
        read_length => {
            is  => 'Boolean',
            doc => 'Whether or not to run read length distribution.',
        },
    ],
    has_metric => [
        _log_directory => {
            is  => 'Text',
            doc => 'Path where workflow logs were written',
        },
        #many other metrics exist--see sub _generate_metrics
    ],
    has => [
        alignment_result => {
            #is => 'Genome::InstrumentData::AlignmentResult::Merged',
            is    => 'Genome::SoftwareResult',
            id_by => 'alignment_result_id',
            doc   => 'the alignment data upon which to run coverage stats',
        },
    ],
    has_transient_optional => [
        log_directory => {
            is  => 'Text',
            doc => 'Path to write logs from running the workflow',
        },
    ],
};

sub resolve_allocation_subdirectory {
    my $self = shift;

    my $hostname = hostname;
    my $user = $ENV{'USER'};
    my $base_dir = sprintf("bam_qc-%s-%s-%s-%s", $hostname, $user, $$, $self->id);

    # TODO: the first subdir is actually specified by the disk management system.
    my $directory = join('/', 'build_merged_alignments','bam_qc',$base_dir);
    return $directory;
}

sub resolve_allocation_disk_group_name {
    $ENV{GENOME_DISK_GROUP_MODELS};
}

sub _staging_disk_usage {
    #need the allocation created in advance for this process
    return 5_000_000; #TODO better estimate
}

sub _working_dir_prefix {
    return 'bam-qc';
}

sub _prepare_staging_directory {
    my $self = shift;

    return $self->temp_staging_directory if ($self->temp_staging_directory);

    unless($self->output_dir) {
        $self->_prepare_output_directory;
    }

    #Stage to network disk because of inner workflow
    my $staging_tempdir = File::Temp->newdir(
        $self->_working_dir_prefix . '-staging-XXXXX',
        DIR     => $self->output_dir,
        CLEANUP => 1,
    );

    $self->temp_staging_directory($staging_tempdir);
}

sub create {
    my $class = shift;
    my $self = $class->SUPER::create(@_);
    return unless ($self);

    $self->_prepare_staging_directory;

    #the bam could be either per lane bam or merged bam, hardcode for now
    my $merge_bam = $self->alignment_result->output_dir.'/'.$self->alignment_result_id . '.bam';
    my $lane_bam  = $self->alignment_result->output_dir.'/'.'all_sequences.bam';

    my ($lane_flag, $bam_file);

    if (-s $merge_bam) {
        $bam_file  = $merge_bam;
    }
    elsif (-s $lane_bam) {
        $bam_file  = $lane_bam;
        $lane_flag = 1;
    }
    else {
        die $self->error_message("Input alignment bam file is missing (no merge bam at $merge_bam or lane bam at $lane_bam)");
    }

    my $fasta_file = $self->alignment_result->reference_build->full_consensus_path('fa');
    die $self->error_message("Reference FASTA File ($fasta_file) is missing") unless -s $fasta_file;
    
    my $log_dir = $self->log_directory;
    unless($log_dir) {
        $log_dir = '' . $self->temp_staging_directory;
    }
    $self->_log_directory($log_dir);

    my %bam_qc_params = (
        output_directory   => '' . $self->temp_staging_directory,
        log_directory      => $log_dir,
        reference_sequence => $fasta_file,
        bam_file           => $bam_file,
        picard_version     => $self->picard_version,
        samtools_version   => $self->samtools_version,
        fastqc_version     => $self->fastqc_version,
        samstat_version    => $self->samstat_version,
        error_rate         => $self->error_rate,
        error_rate_version => $self->error_rate_version,
        read_length        => $self->read_length,
    );

    # use per lane instrument data id as link name
    if ($lane_flag) {
        my $instr_data = $self->alignment_result->instrument_data;
        $bam_qc_params{bam_link_name} = $instr_data->id;
    }

    # Skip samstat html-report for bwamem/sw
    if ($self->alignment_result->aligner_name =~ /^bwa(mem|sw)$/) {
        $bam_qc_params{samstat} = 0;
    }

    my $cmd = Genome::Model::Tools::BamQc::Run->create(%bam_qc_params);
    unless($cmd->execute) {
        die('Failed to run BamQc tool');
    }

    $self->_promote_data;
    $self->_reallocate_disk_allocation;
    $self->_generate_metrics($cmd->output_metrics_hash_ref, $lane_flag);

    return $self;
}

sub _generate_metrics {
    my ($self, $metrics, $flag) = @_;

    $self->status_message('Set bamqc metrics');

    for my $type_label (keys %{$metrics}) {
        # Currently, do not store insert size, quality by cycle, or mean quality histograms
        next if $type_label =~ /Histogram/;
        # Currently, do not store the GcBiasMetrics, 100 Windows of normalized coverage
        next if $type_label eq 'GcBiasMetrics';
        
        my $type_metrics = $metrics->{$type_label};
        # The FlagstatMetrics hashref are one-level
        if ($type_label eq 'FlagstatMetrics') {
            for my $metric_label (keys %{$type_metrics}) {
                my $metric_key = sprintf('bam_qc-%s-%s',$type_label,$metric_label);
                $self->add_metric(metric_name => $metric_key, metric_value => $type_metrics->{$metric_label});
            }
        } 
        else {
            # All other hashrefs are considered to have two-levels, the first level of  
            # the hashref being the key on which lines of metrics are differentiated
            for my $key (keys %{$type_metrics}) {
                for my $metric_label (keys %{$type_metrics->{$key}}) {
                    my $metric_key   = sprintf('bam_qc-%s-%s-%s', $type_label, $key, $metric_label);
                    my $metric_value = $type_metrics->{$key}->{$metric_label};
                    $self->add_metric(metric_name => $metric_key, metric_value => $metric_value);
                }
            }
        }
    }

    # Now fill in alignment metrics
    return 1 unless $flag;
    $self->status_message('Set several alignment metrics');
    
    my %metrics_to_add;
    my $align_result = $self->alignment_result;
    my $instr_data   = $align_result->instrument_data;

    if ($instr_data->is_paired_end) {
        my %convert = (
            read_1_pct_aligned  => ['AlignmentSummaryMetrics', 'CATEGORY-FIRST_OF_PAIR',  'PCT_PF_READS_ALIGNED'],
            read_2_pct_aligned  => ['AlignmentSummaryMetrics', 'CATEGORY-SECOND_OF_PAIR', 'PCT_PF_READS_ALIGNED'],
            read_1_pct_mismatch => ['AlignmentSummaryMetrics', 'CATEGORY-FIRST_OF_PAIR',  'PF_MISMATCH_RATE'],
            read_2_pct_mismatch => ['AlignmentSummaryMetrics', 'CATEGORY-SECOND_OF_PAIR', 'PF_MISMATCH_RATE'],
        );
        %metrics_to_add = $self->_convert_metrics(\%convert, $metrics, 1);
    }
    else { #Only one set of alignment metrics for single_end instrument data
        my %convert = (
            read_1_pct_aligned  => ['AlignmentSummaryMetrics', 'CATEGORY-UNPAIRED', 'PCT_PF_READS_ALIGNED'],
            read_1_pct_mismatch => ['AlignmentSummaryMetrics', 'CATEGORY-UNPAIRED', 'PF_MISMATCH_RATE'],
        );
        %metrics_to_add = $self->_convert_metrics(\%convert, $metrics, 1);
    }
        
    my %convert = (
        median_insert_size => ['InsertSizeMetrics', 'PAIR_ORIENTATION-FR', 'MEDIAN_INSERT_SIZE'],
        sd_insert_size     => ['InsertSizeMetrics', 'PAIR_ORIENTATION-FR', 'STANDARD_DEVIATION'],
    );
    %metrics_to_add = (%metrics_to_add, $self->_convert_metrics(\%convert, $metrics, 0));

    my %align_metrics;
    map{$align_metrics{$_->metric_name} = $_->metric_value}$align_result->metrics;

    for my $metric_key (sort keys %metrics_to_add) {
        if (exists $align_metrics{$metric_key}) {
            $self->warning_message("metric: $metric_key already exist for ".$align_result->id.'. Skip.');
            next;
        }
        $align_result->add_metric(metric_name => $metric_key, metric_value => $metrics_to_add{$metric_key});
    }
    
    return 1;
}

sub _convert_metrics {
    my ($self, $convert, $metrics, $flag) = @_;
    my %convert_metrics;

    for my $metric_name (sort keys %$convert) {
        my $names = $convert->{$metric_name};
        my $metric_value = $metrics->{$names->[0]}->{$names->[1]}->{$names->[2]};
        if ($metric_value) {
            $metric_value = sprintf("%.2f", $metric_value * 100) if $flag;
        }
        else {
            $self->warning_message("Failed to get alignment metric: $metric_name from BamQc output");
            $metric_value = 0;
        }
        $convert_metrics{$metric_name} = $metric_value;
    }

    return %convert_metrics;
}


1;
