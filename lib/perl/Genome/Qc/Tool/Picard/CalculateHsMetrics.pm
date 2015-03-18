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

sub output_file_accessor {
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

1;
