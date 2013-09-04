package Genome::Model::Tools::Analysis::Sammy;

use strict;
use warnings;

use Genome;

class Genome::Model::Tools::Analysis::Sammy {
    is => ['Genome::Model::Tools::Analysis'],
};

sub help_brief {
    "Tools for SAM/BAM analysis with Dan Koboldt's Sammy package.",
}

sub help_synopsis {
    my $self = shift;
    return <<"EOS"
gmt analysis sammy ...
EOS
}

1;

