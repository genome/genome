package Genome::Qc::Tool::Picard::CollectWgsMetrics;

use strict;
use warnings;
use Genome;

class Genome::Qc::Tool::Picard::CollectWgsMetrics {
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
        haploid_coverage => {
            picard_metric => 'MEAN_COVERAGE',
        },
    );
}

sub gmt_class {
    return 'Genome::Model::Tools::Picard::CollectWgsMetrics';
}

1;
