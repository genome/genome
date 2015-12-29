package Genome::Model::MetagenomicShotgun;

use strict;
use warnings;

use Genome;

class Genome::Model::MetagenomicShotgun {
    is => 'Genome::ModelDeprecated',
};

sub do_not_create_define_command { 1 }

sub create {
    die __PACKAGE__ . ' is deprecated.';
}

sub build_subclass_name {
    return 'metagenomic-composition-shotgun';
}

1;

