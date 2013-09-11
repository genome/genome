package Genome::Model::Pair;
use strict;
use warnings;
use Genome;

class Genome::Model::Pair {
    is => 'UR::Value',
    composite_id_separator => '/',
    id_by => [
        first   => { is => 'Genome::Model', id_by => 'first_id' },
        second  => { is => 'Genome::Model', id_by => 'second_id' },
    ],
    doc => 'a pair of genome models',
};

sub __display_name__ {
    return $_[0]->first->__display_name__ . ", " . $_[0]->second->__display_name__;
}

1;

