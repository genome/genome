package Genome::Config::AnalysisMenuItem::Command::List;
use strict;
use warnings;
use Genome;

class Genome::Config::AnalysisMenuItem::Command::List {
    is => 'Genome::Object::Command::List',
    has => [
        subject_class_name  => {
            is_constant => 1, 
            value => 'Genome::Config::AnalysisMenuItem' 
        },
        show => { default_value => 'id,name,configuration_set_id,created_at,updated_at' },
    ],
    doc => 'list analysis menu items',
};

1;

