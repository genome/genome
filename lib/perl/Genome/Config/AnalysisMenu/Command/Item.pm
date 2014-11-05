package Genome::Config::AnalysisMenu::Command::Item;

use strict;
use warnings;

use Genome;

use Genome::Command::Crud;
Genome::Command::Crud->init_sub_commands(
    target_class => 'Genome::Config::AnalysisMenu::Item',
    namespace => 'Genome::Config::AnalysisMenu::Command::Item',
    target_name => 'menu item',
    list => { show => 'id,name,created_by,created_at,updated_at,file_path' },
    update => { include_only => 'description', },
    delete => { do_not_init => 1, },
);

1;
