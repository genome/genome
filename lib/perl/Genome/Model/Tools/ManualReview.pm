package Genome::Model::Tools::ManualReview;

use strict;
use warnings;

use Genome;

class Genome::Model::Tools::ManualReview {
    is => 'Command',
};

sub help_brief {
    "Setup and Run manual review of regions of interest (SNPS, error regions, etc.)"
}

1;
