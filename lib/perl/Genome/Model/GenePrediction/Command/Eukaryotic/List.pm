package Genome::Model::GenePrediction::Command::Eukaryotic::List;

use strict;
use warnings;

use Genome;
use Command; 
use Data::Dumper;

class Genome::Model::GenePrediction::Command::Eukaryotic::List {
    is => 'UR::Object::Command::List',
    has => [
        subject_class_name  => {
            is_constant => 1, 
            value => 'Genome::Model::GenePrediction::Eukaryotic' 
        },
        show => { default_value => 'id,name,subject_name,processing_profile_name,snap_models,fgenesh_model' },
    ],
    doc => 'list genome models',
};

sub sub_command_sort_position { 3 }

1;
