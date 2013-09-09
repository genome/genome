package Genome::Sys::Gateway::Command;

use strict;
use warnings;
use Genome;

class Genome::Sys::Gateway::Command {
    is => 'Command::Tree',
    doc => 'work with other GMS',
};

1;

__END__
Genome::Command::Crud->init_sub_commands(
    target_class => 'Genome::Sys::Service',
    target_name => 'service',
    list => { show => 'name,host,status,pid_status,' },
    create => { do_not_init => 1 },
    update => { do_not_init => 1 },
    delete => { do_not_init => 1 },
);

