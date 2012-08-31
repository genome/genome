package Genome::Model::ImportedAnnotation::Command::ListBuilds;

use strict;
use warnings;

use Genome;
use Command; 
use Data::Dumper;

class Genome::Model::ImportedAnnotation::Command::ListBuilds {
    is => 'UR::Object::Command::List',
    has => [
        subject_class_name  => {
            is_constant => 1, 
            value => 'Genome::Model::Build::ImportedAnnotation' 
        },
        show => { default_value => 'id,name,model.annotation_source,species_name,version,snapshot_date,reference_sequence' },
    ],
    doc => 'list imported annotation builds',
};

sub sub_command_sort_position { 2 }

1;

