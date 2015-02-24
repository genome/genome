package Genome::Interfaces::Comparable;

use strict;
use warnings;
use Genome;

class Genome::Interfaces::Comparable {
    is_abstract => 1,
};

sub compare_output {
    die "Subclass of Genome::Interfaces::Comparable must implement compare_output";
}

1;

