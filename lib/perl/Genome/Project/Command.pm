package Genome::Project::Command;

use strict;
use warnings;

use Genome;

use Genome::Command::Crud;
Genome::Command::Crud->init_sub_commands(
    target_class => 'Genome::Project',
    target_name => 'project',
    list => { show => 'id,name' },
    update => { only_if_null => 1, }
);

1;

