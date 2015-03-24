package Genome::Qc::Tool::PicardCollectMultipleMetrics;

use strict;
use warnings;
use Genome;

class Genome::Qc::Tool::PicardCollectMultipleMetrics {
    is => 'Genome::Qc::Tool',
};

sub supports_streaming {
    return 1;
}

sub cmd_line {
}

sub get_metrics {
}

1;

