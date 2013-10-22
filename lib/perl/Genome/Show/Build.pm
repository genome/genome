package Genome::Show::Build;

use strict;
use warnings;

class Genome::Show::Build {
    is => 'Genome::Show::Base',
};

sub execute {
    my ($self, $object) = @_;
    Genome::Model::Build::Command::View->execute(build => $object);
}
