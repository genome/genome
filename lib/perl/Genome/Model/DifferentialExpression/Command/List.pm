package Genome::Model::DifferentialExpression::Command::List;

use strict;
use warnings;

use Genome;
use Command; 
use Data::Dumper;

class Genome::Model::DifferentialExpression::Command::List {
    is => 'Genome::Object::Command::List',
    has => [
        subject_class_name  => {
            is_constant => 1, 
            value => 'Genome::Model::DifferentialExpression' 
        },
        show => { default_value => 'id,name,processing_profile_id,condition_pairs' },
    ],
    doc => 'list differential expression genome models',
};

1;

