package Genome::Library::Command;

use strict;
use warnings;

use Genome::Command::Crud;
Genome::Command::Crud->init_sub_commands(
    target_class => 'Genome::Library',
    target_name => 'library',
    list => { show => 'id,name,sample_id', },
    update => { only_if_null => 1, exclude => [qw/ sample /], },
    delete => { do_not_init => 1, },
);

1;

