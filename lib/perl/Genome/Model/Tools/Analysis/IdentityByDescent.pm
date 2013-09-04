package Genome::Model::Tools::Analysis::IdentityByDescent;

use strict;
use warnings;

use Genome;

class Genome::Model::Tools::Analysis::IdentityByDescent {
    is => ['Genome::Model::Tools::Analysis'],
};

sub help_brief {
    "Tools to identify identical-by-descent regions in family pedigrees.",
}

sub help_synopsis {
    my $self = shift;
    return <<"EOS"
gmt analysis identity-by-descent --help ...
EOS
}

1;

