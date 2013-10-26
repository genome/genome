package Genome::Model::Tools::Analysis::Illumina;

use strict;
use warnings;

use Genome;

class Genome::Model::Tools::Analysis::Illumina {
    is => ['Genome::Model::Tools::Analysis'],
};

sub help_brief {
    "Tools for analysis of Illumina data.",
}

sub help_synopsis {
    my $self = shift;
    return <<"EOS"
gmt analysis illumina 
EOS
}

1;

