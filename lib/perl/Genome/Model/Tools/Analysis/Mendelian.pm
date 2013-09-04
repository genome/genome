package Genome::Model::Tools::Analysis::Mendelian;

use strict;
use warnings;

use Genome;

class Genome::Model::Tools::Analysis::Mendelian {
    is => ['Genome::Model::Tools::Analysis'],
};

sub help_brief {
    "Tools for analysis of Mendelian disease datasets.",
}

sub help_synopsis {
    my $self = shift;
    return <<"EOS"
gmt analysis mendelian --help ...
EOS
}

1;

