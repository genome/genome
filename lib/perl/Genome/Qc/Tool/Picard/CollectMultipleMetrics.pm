package Genome::Qc::Tool::Picard::CollectMultipleMetrics;

use strict;
use warnings;
use Genome;

class Genome::Qc::Tool::Picard::CollectMultipleMetrics {
    is => 'Genome::Qc::Tool::Picard',
};

sub supports_streaming {
    return 1;
}

sub metrics {
    return (
        alignment_summary_metrics => {
            pct_read_1_aligned => {
                metric_key    => 'CATEGORY-FIRST_OF_PAIR',
                picard_metric => 'PCT_PF_READS_ALIGNED',
            },
            read_1_mismatch_rate => {
                metric_key    => 'CATEGORY-FIRST_OF_PAIR',
                picard_metric => 'PF_MISMATCH_RATE',
            },
            pct_read_2_aligned => {
                metric_key    => 'CATEGORY-SECOND_OF_PAIR',
                picard_metric => 'PCT_PF_READS_ALIGNED',
            },
            read_2_mismatch_rate => {
                metric_key    => 'CATEGORY-SECOND_OF_PAIR',
                picard_metric => 'PF_MISMATCH_RATE',
            },
            pct_properly_paired => {
                metric_key    => 'CATEGORY-PAIR',
                picard_metric => 'PCT_READS_ALIGNED_IN_PAIRS',
            },
            pct_chimeric => {
                metric_key    => 'CATEGORY-PAIR',
                picard_metric => 'PCT_CHIMERAS',
            },
        },
        insert_size_metrics => {
            mean_insert_size => {
                picard_metric => 'MEAN_INSERT_SIZE',
            },
            insert_size_standard_deviation => {
                picard_metric => 'STANDARD_DEVIATION',
            },
            median_insert_size => {
                picard_metric => 'MEDIAN_INSERT_SIZE',
            },
            insert_size_median_absolute_deviation => {
                picard_metric => 'MEDIAN_ABSOLUTE_DEVIATION',
            },
        },
    );
}

sub gmt_class {
    return 'Genome::Model::Tools::Picard::CollectMultipleMetrics';
}

sub gmt_class_for_file_extension {
    my ($self, $file_extension) = @_;
    my %file_extension_to_gmt_class = (
        alignment_summary_metrics    => 'Genome::Model::Tools::Picard::CollectAlignmentSummaryMetrics',
        insert_size_metrics          => 'Genome::Model::Tools::Picard::CollectInsertSizeMetrics',
        quality_by_cycle_metrics     => 'Genome::Model::Tools::Picard::MeanQualityByCycle',
        quality_distribution_metrics => 'Genome::Model::Tools::Picard::QualityScoreDistribution',
    );
    return $file_extension_to_gmt_class{$file_extension};
}

sub qc_metrics_file_accessor {
    return 'output_basename';
}

sub get_metrics {
    my $self = shift;

    my %desired_metric_results;
    my %metrics = $self->metrics;
    while (my ($tool, $tool_metrics) = each %metrics) {
        my $base = $self->qc_metrics_file;
        my $file = join('.', $base, $tool);
        my $gmt_class = $self->gmt_class_for_file_extension($tool);
        my %metric_results = %{$gmt_class->parse_file_into_metrics_hashref($file)};
        my %desired_metric_results_subset = $self->_get_metrics($tool_metrics, \%metric_results);
        @desired_metric_results{ keys %desired_metric_results_subset } = values %desired_metric_results_subset;
    };
    return %desired_metric_results;
}

1;

