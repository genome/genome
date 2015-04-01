package Genome::Qc::Tool::Picard::CollectGcBiasMetrics;

use strict;
use warnings;
use Genome;

class Genome::Qc::Tool::Picard::CollectGcBiasMetrics {
    is => 'Genome::Qc::Tool::Picard',
};

sub supports_streaming {
    return 1;
}

sub output_file_accessor {
    return 'summary_output';
}

sub metrics {
    return (
        at_drop_out => {
            picard_metric => 'AT_DROPOUT',
        },
        gc_drop_out => {
            picard_metric => 'GC_DROPOUT',
        },
    );
}

sub gmt_class {
    return 'Genome::Model::Tools::Picard::CollectGcBiasMetrics';
}

1;
