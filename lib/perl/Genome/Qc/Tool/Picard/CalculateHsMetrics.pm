package Genome::Qc::Tool::Picard::CalculateHsMetrics;

use strict;
use warnings;
use Genome;

class Genome::Qc::Tool::Picard::CalculateHsMetrics {
    is => 'Genome::Qc::Tool::Picard',
};

sub supports_streaming {
    return 1;
}

sub qc_metrics_file_accessor {
    return 'output_file';
}

sub metrics {
    return (
        pct_bases_greater_than_2x_coverage => {
            picard_metric => 'PCT_TARGET_BASES_2X',
        },
        pct_bases_greater_than_10x_coverage => {
            picard_metric => 'PCT_TARGET_BASES_10X',
        },
        pct_bases_greater_than_20x_coverage => {
            picard_metric => 'PCT_TARGET_BASES_20X',
        },
        pct_bases_greater_than_30x_coverage => {
            picard_metric => 'PCT_TARGET_BASES_30X',
        },
        pct_bases_greater_than_40x_coverage => {
            picard_metric => 'PCT_TARGET_BASES_40X',
        },
        pct_bases_greater_than_50x_coverage => {
            picard_metric => 'PCT_TARGET_BASES_50X',
        },
        pct_bases_greater_than_100x_coverage => {
            picard_metric => 'PCT_TARGET_BASES_100X',
        },
    );
}

sub gmt_class {
    return 'Genome::Model::Tools::Picard::CalculateHsMetrics';
}

sub bait_intervals {
    my $self = shift;

    my $translated_list = $self->translate_feature_list('tiled_region');
    return $self->create_interval_file($translated_list);
}

sub target_intervals {
    my $self = shift;

    my $translated_list = $self->translate_feature_list('target_region');
    return $self->create_interval_file($translated_list);
}

sub translate_feature_list {
    my ($self, $track_name) = @_;

    my $translated_list = Genome::Sys->create_temp_file_path;
    my %params = (
        feature_list => $self->alignment_result->instrument_data->target_region_set,
        track_name => $track_name,
        output_path => $translated_list,
        merge => 0,
    );
    unless ($self->reference_build eq $self->alignment_result->instrument_data->target_region_set->reference) {
        $params{alternate_reference} = $self->reference_build;
    }
    my $translated_list_cmd = Genome::FeatureList::Command::DumpMergedList->create(%params);
    $translated_list_cmd->execute;

    return $translated_list;
}

sub create_interval_file {
    my ($self, $translated_list) = @_;

    my $interval_file = Genome::Sys->create_temp_file_path;
    my $interval_file_cmd = Genome::Model::Tools::Bed::ToIntervals->create(
        bed_file => $translated_list,
        interval_file => $interval_file,
        seqdict_file => $self->reference_build->get_sequence_dictionary('sam'),
    );
    $interval_file_cmd->execute;

    return $interval_file;
}

1;
