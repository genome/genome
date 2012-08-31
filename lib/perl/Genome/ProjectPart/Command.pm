package Genome::ProjectPart::Command;

use strict;
use warnings;

use Genome;

class Genome::ProjectPart::Command {
    is => 'Command::Tree',
    doc => 'work with project parts',
};

Genome::Command::Crud->init_sub_commands(
    target_class => 'Genome::ProjectPart',
    target_name => 'project-part',
    list => { show => 'project_id,entity_class_name,entity_id,label,role' },
    update => { only_if_null => 1, },
);

1;


