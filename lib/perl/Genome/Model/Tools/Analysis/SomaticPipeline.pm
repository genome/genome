package Genome::Model::Tools::Analysis::SomaticPipeline;

use strict;
use warnings;

use Genome;

class Genome::Model::Tools::Analysis::SomaticPipeline {
    is => ['Genome::Model::Tools::Analysis'],
};

sub help_brief {
    "Tools related to the genome model somatic pipeline.",
}

sub help_synopsis {
    my $self = shift;
    return <<"EOS"
gmt analysis somatic-pipeline ...
EOS
}

1;

