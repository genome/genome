package Genome::InstrumentData::AlignmentResult::Merged::Tophat2AlignmentStats;

use strict;
use warnings;

use Genome;
use Sys::Hostname;
use File::Path;

class Genome::InstrumentData::AlignmentResult::Merged::Tophat2AlignmentStats {
    is => ['Genome::SoftwareResult::Stageable'],
    has_input => [
        alignment_result_id => {
            is => 'Number',
            doc => 'ID of the result for the alignment data upon which to run coverage stats',
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
    my $base_dir = sprintf("tophat2_alignment_stats-%s-%s-%s-%s", $hostname, $user, $$, $self->id);

    # TODO: the first subdir is actually specified by the disk management system.
    my $directory = join('/', 'build_merged_alignments','tophat2_alignment_stats',$base_dir);
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
    return 'tophat2-alignment-stats';
}

sub create {
    my $class = shift;
    my $self = $class->SUPER::create(@_);
    return unless ($self);

    $self->_prepare_staging_directory;

    my $bam_file = $self->alignment_result->bam_file;
    my $alignment_stats_file = $self->temp_staging_directory .'/alignment_stats.txt';

    unless (Genome::Model::Tools::BioSamtools::Tophat2AlignmentStats->execute(
        bam_file => $bam_file,
        alignment_stats_file => $alignment_stats_file,
    )) {
        $self->error_message('Failed to run BioSamtools Tophat2AlignmentStats for bam file: '. $bam_file);
        return;
    }
    unless ($self->output_dir) {
        $self->_prepare_output_directory;
    }
    $self->_promote_data;
    $self->_reallocate_disk_allocation;

    $self->_generate_metrics();

    return $self;
}

sub _generate_metrics {
    my $self = shift;

    my $alignment_stats_file = $self->output_dir .'/alignment_stats.txt';
    my $metrics_hash_ref = Genome::Model::Tools::BioSamtools::Tophat2AlignmentStats->parse_alignment_stats_summary_hash_ref($alignment_stats_file);
    for my $metric_name (keys %{$metrics_hash_ref}) {
        $self->add_metric(
            metric_name => $metric_name,
            metric_value => $metrics_hash_ref->{$metric_name},
        );
    }
    
    return 1;
}


1;
