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

sub qc_metrics_file_accessor {
    return 'summary_output';
}

sub gmt_class {
    return 'Genome::Model::Tools::Picard::CollectGcBiasMetrics';
}

sub chart_output {
    return Genome::Sys->create_temp_file_path();
}

sub output_file {
    return Genome::Sys->create_temp_file_path();
}

1;
