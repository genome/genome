package Genome::Model::MutationalSignificance::Command::List;

use strict;
use warnings;

use Genome;
use Command; 
use Data::Dumper;

class Genome::Model::MutationalSignificance::Command::List {
    is => 'Genome::Object::Command::List',
    has => [
        subject_class_name  => {
            is_constant => 1, 
            value => 'Genome::Model::MutationalSignificance' 
        },
        show => { default_value => 'id,name' },
    ],
    doc => 'list mutational-significance genome models',
};

1;

