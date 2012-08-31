package Genome::Individual::Command;

use strict;
use warnings;

use Genome::Command::Crud;
Genome::Command::Crud->init_sub_commands(
    target_class => 'Genome::Individual',
    target_name => 'individual',
    list => { show => 'id,name,upn,species_name,common_name,gender', },
    update => { only_if_null => 1, },
    delete => { do_not_init => 1, },
);

1;

