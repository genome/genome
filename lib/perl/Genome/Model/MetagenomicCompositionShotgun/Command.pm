package Genome::Model::MetagenomicCompositionShotgun::Command;

use Genome;
use strict;
use warnings;


class Genome::Model::MetagenomicCompositionShotgun::Command {
    is => 'Command',
    doc => 'operate on metagenomic-composition-shotgun models',
};

sub sub_command_category { 'type specific' }

1;

