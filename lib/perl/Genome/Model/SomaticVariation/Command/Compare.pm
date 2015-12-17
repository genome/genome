use Genome;
use strict;
use warnings;

package Genome::Model::SomaticVariation::Command::Compare;

class Genome::Model::SomaticVariation::Command::Compare {
    is => 'Command::Tree',
    doc => 'compare a pair of somatic variation builds in a variety of ways',
};

1;

