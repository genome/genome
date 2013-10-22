package Genome::Show::Base;

use strict;
use warnings;

class Genome::Show::Base {
    is_abstract => 1,
};

# Must be implemented
sub execute {
    my ($self, $object) = @_;
}
