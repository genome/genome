package Genome::Qc::Command::Config::List;
use strict;
use warnings;
use Genome;

class Genome::Qc::Command::Config::List {
    is => 'Genome::Object::Command::List',
    has => [
        subject_class_name  => {
            is_constant => 1,
            value => 'Genome::Qc::Config'
        },
        show => { default_value => 'id,name,type' },
    ],
    doc => 'List Qc configurations',
};

sub sub_command_sort_position { 3 }

1;

