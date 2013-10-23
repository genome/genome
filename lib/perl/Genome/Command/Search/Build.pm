package Genome::Command::Search::Build;

use strict;
use warnings;

class Genome::Command::Search::Build {
    is => 'Genome::Command::Search::Base',
};

sub display_single {
    my ($self, $object) = @_;
    Genome::Model::Build::Command::View->execute(build => $object);
}

sub display_many {
    my $self = shift;
    Genome::Model::Build::Command::Status->execute(builds => \@_);
}
