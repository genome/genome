package Genome::Model::Metric::Command::List;

use strict;
use warnings;

use Genome;

class Genome::Model::Metric::Command::List {
    is => 'Genome::Model::Command::BuildRelatedList',
    has => [
        subject_class_name  => {
            is_constant => 1, 
            value => 'Genome::Model::Metric' 
        },
        show => { default_value => 'model_id,build_id,name,value' },
    ],
    doc => 'list genome model per-build metrics',
};

sub sub_command_sort_position { 1 }

1;

