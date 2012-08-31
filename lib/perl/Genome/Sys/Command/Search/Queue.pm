package Genome::Sys::Command::Search::Queue;

use strict;
use warnings;

use Genome;

class Genome::Sys::Command::Search::Queue {
    is => 'Command::Tree',
    doc => 'work with search queue',
};

use Genome::Command::Crud;
Genome::Command::Crud->init_sub_commands(
    target_class => 'Genome::Search::Queue',
    namespace => 'Genome::Sys::Command::Search::Queue',
    target_name => 'queue',
    list => {
        show => 'priority,timestamp,subject_class,subject_id',
        order_by => 'priority,timestamp',
    },
);

1;
