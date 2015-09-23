package Genome::Qc::Tool::Picard::CollectAlignmentSummaryMetrics;

use strict;
use warnings;
use Genome;

class Genome::Qc::Tool::Picard::CollectAlignmentSummaryMetrics {
    is => 'Genome::Qc::Tool::Picard',
};

sub supports_streaming {
    return 1;
}

sub qc_metrics_file_accessor {
    return 'output_file';
}

sub gmt_class {
    return 'Genome::Model::Tools::Picard::CollectAlignmentSummaryMetrics';
}

1;
