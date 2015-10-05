package Genome::Qc::Tool::Picard::CollectInsertSizeMetrics;

use strict;
use warnings;
use Genome;

class Genome::Qc::Tool::Picard::CollectInsertSizeMetrics {
    is => 'Genome::Qc::Tool::Picard',
};

sub supports_streaming {
    return 1;
}

sub qc_metrics_file_accessor {
    return 'output_file';
}

sub gmt_class {
    return 'Genome::Model::Tools::Picard::CollectInsertSizeMetrics';
}

sub histogram_file {
    return Genome::Sys->create_temp_file_path;
}

1;
