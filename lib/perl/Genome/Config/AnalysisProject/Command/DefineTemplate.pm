package Genome::Config::AnalysisProject::Command::DefineTemplate;

use strict;
use warnings;

class Genome::Config::AnalysisProject::Command::DefineTemplate {
    is => 'Genome::Config::AnalysisProject::Command::Create',
};

sub help_brief {
    return 'Define a template for an analysis-project';
}

sub help_synopsis {
    return 'genome config analysis-project define-template "my-template"';
}

sub help_detail {
    return <<EOS;
Setup a new "template" analysis project.  This will not process
instrument data, but can be used to store a collection of configuration
to be copied later.
EOS
}

sub execute {
    my $self = shift;

    my $project = $self->SUPER::_execute_body;
    return unless $project;

    $project->status('Template');
    $self->status_message('Converted new project to a Template. It cannot be used to create models.');
    return $project;
}

1;
