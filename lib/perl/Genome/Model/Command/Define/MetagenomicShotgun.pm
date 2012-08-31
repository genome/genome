package Genome::Model::Command::Define::MetagenomicShotgun;

use strict;
use warnings;

use Genome;

class Genome::Model::Command::Define::MetagenomicShotgun {
    is => 'Genome::Model::Command::Define::Base',
    has => [
        _model_class => { value => 'Genome::Model::MetagenomicShotgun', },
    ],
    doc => 'define a new metagenomic shotgun genome model',
};


1; 
