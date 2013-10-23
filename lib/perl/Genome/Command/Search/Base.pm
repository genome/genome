package Genome::Command::Search::Base;

use strict;
use warnings;

class Genome::Command::Search::Base {
    is_abstract => 1,
};

# Must be implemented
sub display_single {
    my ($self, $object) = @_;
}

sub display_many {
    my $self = shift;
    my @objects = @_;
}
