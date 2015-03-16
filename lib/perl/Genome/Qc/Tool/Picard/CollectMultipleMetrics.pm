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

sub cmd_line {
}

1;

