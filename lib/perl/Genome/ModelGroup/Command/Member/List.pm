package Genome::ModelGroup::Command::Member::List;

use strict;
use warnings;

use Genome;

class Genome::ModelGroup::Command::Member::List {
    is => 'UR::Object::Command::List',
    has => [
        show => {
            doc => 'properties of the member models to list (comma-delimited)',
            is_optional => 1,
            default_value => 'model.id,model.name',
        },
        subject_class_name  => {
            is_constant => 1,
            value => 'Genome::ModelGroupBridge'
        },
    ],
    doc => 'list the member models of a model-group',
};

1;

