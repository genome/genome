package Genome::Command::Search::Build;

use strict;
use warnings;

class Genome::Command::Search::Build {
    is => 'Genome::Command::Search::Base',
};

sub execute {
    my ($self, $object) = @_;
    Genome::Model::Build::Command::View->execute(build => $object);
}
