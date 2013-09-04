package Genome::Config::AnalysisProject::SubjectPairing::Command;

use strict;
use warnings;

use Genome::Command::Crud;
Genome::Command::Crud->init_sub_commands(
    target_class => 'Genome::Config::AnalysisProject::SubjectPairing',
    target_name => 'subject_pairing',
    list => { show => 'id,analysis_project,control_subject,experimental_subject,created_by,created_at,updated_at' },
    create => { exclude => ['created_at', 'updated_at'] },
    update => { do_not_init => 1 },
);

class Genome::Config::AnalysisProject::SubjectPairing::Command {
    is => 'Command::Tree',
    doc => 'commands that deal with subject pairing for analysis projects',
};

1;
