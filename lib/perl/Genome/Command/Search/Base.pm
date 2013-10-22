package Genome::Command::Search::Base;

use strict;
use warnings;

class Genome::Command::Search::Base {
    is_abstract => 1,
};

# Must be implemented
sub execute {
    my ($self, $object) = @_;
}
