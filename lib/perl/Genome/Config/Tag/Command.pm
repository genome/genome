package Genome::Config::Tag::Command;

use strict;
use warnings;

use Genome::Command::Crud;
Genome::Command::Crud->init_sub_commands(
    target_class => 'Genome::Config::Tag',
    target_name => 'tag',
    list => { show => 'id,name,description' },
    create => { do_not_init => 1 },
    update => { include_only => ['description'] },
    delete => { do_not_init => 1 },
);

class Genome::Config::Tag::Command {
    is => 'Command::Tree',
    doc => 'commands to manipulate configuration tags',
};

1;
