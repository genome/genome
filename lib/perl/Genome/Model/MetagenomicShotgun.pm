package Genome::Model::MetagenomicShotgun;

use strict;
use warnings;

use Genome;

class Genome::Model::MetagenomicShotgun {
    is => 'Genome::ModelDeprecated',
};

sub create {
    die __PACKAGE__ . ' is deprecated.';
}

sub build_subclass_name {
    return 'metagenomic-composition-shotgun';
}

1;

