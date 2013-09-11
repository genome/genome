package Genome::Model::Pair;
use strict;
use warnings;
use Genome;

class Genome::Model::Pair {
    is => 'UR::Value',
    id_by => [
        first   => { is => 'Genome::Model', id_by => 'first_id' },
        second  => { is => 'Genome::Model', id_by => 'second_id' },
    ],
    doc => 'a pair of genome models',
};

1;

