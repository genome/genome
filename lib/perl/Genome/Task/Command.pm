package Genome::Task::Command;

use strict;
use warnings;

use Genome;

class Genome::Task::Command {
    is => 'Command::Tree',
    doc => 'work with tasks',
};

Genome::Command::Crud->init_sub_commands(
    target_class => 'Genome::Task',
    target_name => 'task',
    list => { show => 'id,command_class,status,user_id,time_submitted,time_started,time_finished' },
    create => { do_not_init => 1 },
    update => { do_not_init => 1 },
    delete => { do_not_init => 1 },
);


1;
