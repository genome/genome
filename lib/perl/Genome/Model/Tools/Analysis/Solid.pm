package Genome::Model::Tools::Analysis::Solid;

use strict;
use warnings;

use Genome;

class Genome::Model::Tools::Analysis::Solid {
    is => ['Genome::Model::Tools::Analysis'],
};

sub help_brief {
    "Tools for analysis of ABI/Solid data.",
}

sub help_synopsis {
    my $self = shift;
    return <<"EOS"
genome-model tools analysis solid ...
EOS
}

1;

