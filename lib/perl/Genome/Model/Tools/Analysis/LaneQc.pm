package Genome::Model::Tools::Analysis::LaneQc;

use strict;
use warnings;

use Genome;

class Genome::Model::Tools::Analysis::LaneQc {
    is => ['Genome::Model::Tools::Analysis'],
};

sub help_brief {
    "Tools for QC-checks of Illumina/Solexa data.",
}

sub help_synopsis {
    my $self = shift;
    return <<"EOS"
gmt analysis lane-qc --help ...
EOS
}

1;

