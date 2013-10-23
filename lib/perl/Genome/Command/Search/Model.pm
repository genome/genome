package Genome::Command::Search::Model;

use strict;
use warnings;

class Genome::Command::Search::Model {
    is => 'Genome::Command::Search::Base',
};

sub display_single {
    my ($self, $object) = @_;
    my @builds = $object->builds;
    Genome::Model::Build::Command::Status->execute(builds => \@builds);
}

sub display_many {
    my $self = shift;
    Genome::Model::Command::Status->execute(models => \@_);
}
