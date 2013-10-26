package Genome::Model::Tools::Analysis::CaseControl;

use strict;
use warnings;

use Genome;

class Genome::Model::Tools::Analysis::CaseControl {
    is => ['Genome::Model::Tools::Analysis'],
};

sub help_brief {
    "Tools for analysis of CaseControl disease datasets.",
}

sub help_synopsis {
    my $self = shift;
    return <<"EOS"
gmt analysis case-control --help ...
EOS
}

1;

