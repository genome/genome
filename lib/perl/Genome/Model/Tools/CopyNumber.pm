package Genome::Model::Tools::CopyNumber;

use strict;
use warnings;

use Genome;

class Genome::Model::Tools::CopyNumber {
    is => ['Command'],
};

sub help_brief {
    "Tools and R scripts to perform copy number analysis."
}

sub help_detail {
    return <<EOS
Tools and R scripts to perform copy number analysis.
EOS
}

1;
