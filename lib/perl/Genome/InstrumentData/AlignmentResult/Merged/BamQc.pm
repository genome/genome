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
            is => 'Number',
            doc => 'ID of the result for the alignment data upon which to run coverage stats',
        },
    ],
    has_param => [
        picard_version => {
            is => 'Text',
            doc => 'The version of Picard to use.',
        },
        samstat_version => {
            is => 'Text',
            doc => 'The version of SamStat to use.',
        },
        fastqc_version => {
            is => 'Text',
            doc => 'The version of FastQC to use.',
        },
        samtools_version => {
            is => 'Text',
            doc => 'The version of Samtools to use.',
        },
        error_rate => {
            is => 'Boolean',
            doc => 'Whether or not to generate error rate summaries.',
        },
        error_rate_pileup => {
            is => 'Boolean',
            doc => 'Whether or not to iterate over every position to calculate positional error rates.',
        },
    ],
    has_metric => [
        _log_directory => {
            is => 'Text',
            doc => 'Path where workflow logs were written',
        },
        #many other metrics exist--see sub _generate_metrics
    ],
    has => [
        alignment_result => {
            is => 'Genome::InstrumentData::AlignmentResult::Merged',
            id_by => 'alignment_result_id',
            doc => 'the alignment data upon which to run coverage stats',
        },
    ],
    has_transient_optional => [
        log_directory => {
            is => 'Text',
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
    return 'info_genome_models';
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

sub _gather_params_for_get_or_create {
    my $class = shift;

    my $bx = UR::BoolExpr->resolve_normalized_rule_for_class_and_params($class, @_);

    my %params = $bx->params_list;
    my %is_input;
    my %is_param;
    my $class_object = $class->__meta__;
    for my $key ($class->property_names) {
        my $meta = $class_object->property_meta_for_name($key);
        if ($meta->{is_input} && exists $params{$key}) {
            $is_input{$key} = $params{$key};
        } elsif ($meta->{is_param} && exists $params{$key}) {
            $is_param{$key} = $params{$key};
        }
    }

    my $inputs_bx = UR::BoolExpr->resolve_normalized_rule_for_class_and_params($class, %is_input);
    my $params_bx = UR::BoolExpr->resolve_normalized_rule_for_class_and_params($class, %is_param);

    my %software_result_params = (
        params_id=>$params_bx->id,
        inputs_id=>$inputs_bx->id,
        subclass_name=>$class
    );

    return {
        software_result_params => \%software_result_params,
        subclass => $class,
        inputs=>\%is_input,
        params=>\%is_param,
    };
}

sub create {
    my $class = shift;
    my $self = $class->SUPER::create(@_);
    return unless ($self);

    $self->_prepare_staging_directory;

    my $fasta_file = $self->alignment_result->reference_build->full_consensus_path('fa');
    my $bam_file = $self->alignment_result->merged_alignment_bam_path;

    die $self->error_message("Reference FASTA File ($fasta_file) is missing") unless -s $fasta_file;
    die $self->error_message("Bam File ($bam_file) is missing") unless -s $bam_file;

    my $log_dir = $self->log_directory;
    unless($log_dir) {
        $log_dir = '' . $self->temp_staging_directory;
    }
    $self->_log_directory($log_dir);

    my %bam_qc_params = (
        output_directory => '' . $self->temp_staging_directory,
        log_directory => $log_dir,
        reference_sequence => $fasta_file,
        bam_file => $bam_file,
        picard_version => $self->picard_version,
        samtools_version => $self->samtools_version,
        fastqc_version => $self->fastqc_version,
        samstat_version => $self->samstat_version,
        error_rate =>  $self->error_rate,
        error_rate_pileup => $self->error_rate_pileup,
    );

    my $cmd = Genome::Model::Tools::BamQc::Run->create(%bam_qc_params);
    unless($cmd->execute) {
        die('Failed to run BamQc tool');
    }

    $self->_promote_data;
    $self->_reallocate_disk_allocation;

    $self->_generate_metrics($cmd->output_metrics_hash_ref);

    return $self;
}

sub _generate_metrics {
    my $self = shift;
    my $metrics = shift;
    
    for my $type_label (keys %{$metrics}) {
        # Currently, do not store insert size, quality by cycle, or mean quality histograms
        if ($type_label =~ /Histogram/) { next; }
        # Currently, do not store the GcBiasMetrics, 100 Windows of normalized coverage
        if ($type_label eq 'GcBiasMetrics') { next; }
        
        my $type_metrics = $metrics->{$type_label};
        # The FlagstatMetrics hashref are one-level
        if ($type_label eq 'FlagstatMetrics') {
            for my $metric_label (keys %{$type_metrics}) {
                my $metric_key = sprintf('bam_qc-%s-%s',$type_label,$metric_label);
                $self->add_metric(metric_name => $metric_key, metric_value => $type_metrics->{$metric_label});
            }
        } else {
            # All other hashrefs are considered to habe two-levels, the first level of the hashref being the key on which lines of metrics are differentiated
            for my $key (keys %{$type_metrics}) {
                for my $metric_label (keys %{$type_metrics->{$key}}) {
                    my $metric_key = sprintf('bam_qc-%s-%s-%s',$type_label,$key,$metric_label);
                    $self->add_metric(metric_name => $metric_key, metric_value => $type_metrics->{$key}->{$metric_label});
                }
            }
        }
    }

    return 1;
}


1;
