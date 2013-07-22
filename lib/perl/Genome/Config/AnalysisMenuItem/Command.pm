package Genome::Config::AnalysisMenuItem::Command;

use strict;
use warnings;

use Genome::Command::Crud;
Genome::Command::Crud->init_sub_commands(
    target_class => 'Genome::Config::AnalysisMenuItem',
    target_name => 'analysis_project',
    list => { show => 'id,name,configuration_set_id,created_at,updated_at' },
    create => { do_not_init => 1 },
    update => { do_not_init => 1 },

);

class Genome::Config::AnalysisMenuItem::Command {
    is => 'Command::Tree',
    doc => 'commands that deal with analysis menu items',
};

1;
