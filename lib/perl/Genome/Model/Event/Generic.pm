package Genome::Model::Event::Generic;

use strict;
use warnings;

use Genome;

class Genome::Model::Event::Generic {
    is => ['Genome::Model::Event'],
};

sub help_brief {
    "Generic placeholder class for old events that have been removed from the code base."
}

sub help_detail {
    return <<"EOS"
    Generic placeholder class for old events that have been removed from the code base.
EOS
}

1;
