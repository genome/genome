package Genome::Config::AnalysisMenu::Item::Command::List;

use strict;
use warnings;
use Genome;

class Genome::Config::AnalysisMenu::Item::Command::List {
    is => 'Genome::Object::Command::List',
    has => [
        subject_class_name => {
            is_constant => 1,
            value => 'Genome::Config::AnalysisMenu::Item'
        },
        show => {default_value => 'id,name,created_by,created_at,updated_at,file_path'},
    ],
    doc => 'list Analysis Menu Items',
};

sub sub_command_sort_position { 1 }

1;
