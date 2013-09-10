package Genome::Model::Pair;
use strict;
use warnings;
use Genome;

class Genome::Model::Pair {
    is => 'UR::Value',
    id_by => [
        first   => { is => 'Genome::Model', id_by => 'id1' },
        second  => { is => 'Genome::Model', id_by => 'id2' },
    ],
    doc => 'a pair of genome models',
};

1;

