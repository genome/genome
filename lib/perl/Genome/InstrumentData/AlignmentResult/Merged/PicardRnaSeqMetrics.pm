package Genome::InstrumentData::AlignmentResult::Merged::PicardRnaSeqMetrics;

use strict;
use warnings;

use Genome;
use Sys::Hostname;
use File::Path;

class Genome::InstrumentData::AlignmentResult::Merged::PicardRnaSeqMetrics {
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
            doc => 'The version of cufflinks to use.',
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
    my $base_dir = sprintf("picard_rna_seq_metrics-%s-%s-%s-%s", $hostname, $user, $$, $self->id);

    # TODO: the first subdir is actually specified by the disk management system.
    my $directory = join('/', 'build_merged_alignments','picard_rna_seq_metrics',$base_dir);
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
    return 'picard-rna-seq-metrics';
}


sub create {
    my $class = shift;
    my $self = $class->SUPER::create(@_);
    return unless ($self);

    $self->_prepare_staging_directory;

    my $alignment_result = $self->alignment_result;

    my $log_dir = $self->log_directory;
    unless($log_dir) {
        $log_dir = '' . $self->temp_staging_directory;
    }
    $self->_log_directory($log_dir);

    my $metrics_directory = $self->temp_staging_directory;
    unless (Genome::InstrumentData::AlignmentResult::Command::PicardRnaSeqMetrics->execute(
        metrics_directory => $metrics_directory,
        annotation_build => $alignment_result->annotation_build,
        reference_build => $alignment_result->reference_build,
        alignment_result => $alignment_result,
        picard_version => $self->picard_version,
    )) {
        return;
    }
    
    $self->_prepare_output_directory;
    $self->_promote_data;
    $self->_reallocate_disk_allocation;

    $self->_generate_metrics();

    return $self;
}

sub _generate_metrics {
    my $self = shift;

    my $rna_seq_metrics = $self->output_dir .'/'. Genome::Model::Build::RnaSeq->default_metrics_output_filename;
    my $metrics_hash_ref = Genome::Model::Tools::Picard::CollectRnaSeqMetrics->parse_file_into_metrics_hashref($rna_seq_metrics);
    for my $metric_name (keys %{$metrics_hash_ref}) {
        $self->add_metric(metric_name => $metric_name, metric_value => $metrics_hash_ref->{$metric_name});
    }
    return 1;
}


1;
