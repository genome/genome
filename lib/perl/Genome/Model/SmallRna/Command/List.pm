package Genome::Model::SmallRna::Command::List;

use strict;
use warnings;

use Genome;
use Command;
use Data::Dumper;

class Genome::Model::SmallRna::Command::List {
    is => 'UR::Object::Command::List',
    has => [
        subject_class_name  => {
            is_constant => 1,
            value => 'Genome::Model::SmallRna'
        },
        show => { default_value => 'id,name,ref_model' },
    ],
    doc => 'list small-rna genome models',
};

1;

