package Genome::Config::AnalysisProject::Command::List;
use strict;
use warnings;
use Genome;

class Genome::Config::AnalysisProject::Command::List {
    is => 'Genome::Object::Command::List',
    has => [
        subject_class_name  => {
            is_constant => 1, 
            value => 'Genome::Config::AnalysisProject'
        },
        show => { default_value => 'id,name,created_by,_configuration_set_id,created_at,updated_at' },
    ],
    doc => 'list analysis projects',
};

1;

