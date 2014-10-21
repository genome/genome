package Genome::Config::AnalysisProject::Command::AddModelToAnalysisProject;

use strict;
use warnings;

use Genome;

class Genome::Config::AnalysisProject::Command::AddModelToAnalysisProject {
    is => 'Command::V2',
    has_input => [
       analysis_project => {
            is                  => 'Genome::Config::AnalysisProject',
            doc                 => 'the analysis project to which to add the model',
            shell_args_position => 1,
       },
       models  => {
            is                  => 'Genome::Model',
            doc                 => 'model(s) to add to this analysis project',
            is_many             => 1,
            shell_args_position => 2,
       },
    ],
    has_optional_input => {
        allow_projects_not_in_progress => {
            is                  => 'Boolean',
            doc                 => 'add the model(s) even if the analysis project is not "In Progress"',
            default_value       => 0,
        },
    },
};

sub help_brief {
    return 'associate a model with an analysis project';
}

sub help_synopsis {
    return 'genome config analysis-project add-model-to-analysis-project <analysis-project> <model> [<model> ...]';
}

sub help_detail {
    return <<EOS
This will associate one or more models with an analysis project.  This will fail if a model is already associated with an analysis project.
EOS
}

sub execute {
    my $self = shift;

    my $analysis_project = $self->analysis_project;
    return unless $self->_verify_analysis_project_status($analysis_project);

    my @models = $self->models;
    return unless $self->_verify_models_are_assignable(@models);

    for my $model (@models) {
        $self->_assign_model($analysis_project, $model);
    }

    $self->status_message('Models added to analysis project.');

    return 1;
}

sub _verify_analysis_project_status {
    my $self = shift;
    my $analysis_project = shift;

    unless($analysis_project->status eq 'In Progress') {
        if($self->allow_projects_not_in_progress) {
            $self->debug_message(
                'Allowing assignment to Analysis Project %s in status %s',
                $analysis_project->id,
                $analysis_project->status
            );
        } else {
            $self->error_message(
                'Analysis Project %s has status %s. Not assigning models.',
                $analysis_project->id,
                $analysis_project->status
            );
            return;
        }
    }

    return 1;
}

sub _verify_models_are_assignable {
    my $self = shift;
    my @models = @_;

    my $analysis_project = $self->analysis_project;

    my $all_models_are_assignable = 1;
    for my $model (@models) {
        if($model->analysis_projects) {
            my $model_anp = $model->analysis_projects;
            unless($model_anp eq $analysis_project) {
                $self->error_message(
                    'Model %s is already assigned to analysis project %s.  Cannot assign it to %s.',
                    $model->id,
                    $model_anp->id,
                    $analysis_project->id,
                );
                $all_models_are_assignable = 0;
            }
        }
    }

    return $all_models_are_assignable;
}

sub _assign_model {
    my $self = shift;
    my $analysis_project = shift;
    my $model = shift;

    if($model->analysis_projects) {
        $self->debug_message(
            'Model %s is already assigned to analysis project %s.',
            $model->id,
            $model->analysis_projects->id,
        );
        return 1;
    }

    $analysis_project->add_model_bridge(model => $model);
    $self->debug_message(
        'Model %s assigned to analysis project %s.',
        $model->id,
        $model->analysis_projects->id,
    );

    return 1;
}

1;
