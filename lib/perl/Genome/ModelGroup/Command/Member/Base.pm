package Genome::ModelGroup::Command::Member::Base;

use strict;
use warnings;

use Genome;

class Genome::ModelGroup::Command::Member::Base {
    is => 'Command::V2',
    is_abstract => 1,
    has => [
        model_group => { 
            is => 'Genome::ModelGroup', 
            shell_args_position => 1,
            doc => 'Model group name or id.',
        },
    ],
    doc => "work with the members of model-groups",
};

1;

