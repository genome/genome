package Genome::Model::Tools::Convergence;

use strict;
use warnings;

use Genome;

class Genome::Model::Tools::Convergence {
    is => 'Command',
    has => [
    ],
};

sub help_brief {
    "Tools to work with harmonic convergence.",
}

sub help_synopsis {
    my $self = shift;
    return <<"EOS"
gmt convergence ...
EOS
}

sub help_detail {                           
    return <<EOS 
Tools to run Harmonic Convergence analysis on a set of models.
EOS
}

1;
