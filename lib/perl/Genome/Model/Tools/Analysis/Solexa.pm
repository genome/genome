package Genome::Model::Tools::Analysis::Solexa;

use strict;
use warnings;

use Genome;

class Genome::Model::Tools::Analysis::Solexa {
    is => ['Genome::Model::Tools::Analysis'],
};

sub help_brief {
    "DEPRECATED tools for analysis of Illumina/Solexa data.",
}

sub help_synopsis {
    my $self = shift;
    return <<"EOS"
gmt analysis solexa ...
EOS
}

1;

