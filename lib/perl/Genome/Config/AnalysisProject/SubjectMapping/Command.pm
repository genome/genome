package Genome::Config::AnalysisProject::SubjectMapping::Command;

use strict;
use warnings;

use Genome::Command::Crud;

class Genome::Config::AnalysisProject::SubjectMapping::Command {
    is => 'Command::Tree',
    doc => 'commands that deal with subject mapping for analysis projects',
};

Genome::Command::Crud->init_sub_commands(
    target_class => 'Genome::Config::AnalysisProject::SubjectMapping',
    target_name => 'subject-mapping',
    list => { show => 'id,subjects,tags.name' },
    update => { do_not_init => 1, },
    create => { do_not_init => 1, },
);

1;
