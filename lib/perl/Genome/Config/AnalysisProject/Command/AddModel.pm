package Genome::Config::AnalysisProject::Command::AddModel;

use strict;
use warnings;

use Genome;

class Genome::Config::AnalysisProject::Command::AddModel {
    is => 'Command::V2',
    has_input => [
       profile_item => {
            is                  => 'Genome::Config::Profile::Item',
            doc                 => 'the configuration to which to associate the model(s) to the analysis project',
            shell_args_position => 1,
       },
       models  => {
            is                  => 'Genome::Model',
            doc                 => 'model(s) to add to the analysis project',
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
    return 'genome config analysis-project add-model <profile-item> <model> [<model> ...]';
}

sub help_detail {
    return <<EOS
This will associate one or more models with an analysis project.  This will fail if a model is already associated with an analysis project.

This command expects a manual configuration to already exist for the model(s) to be added.  The available configurations for an analysis-project can be listed with `genome analysis-project show-config`.  If one does not already exist, it can be added using `genome analysis-project add-config-file` with the --store-only option.
EOS
}

sub execute {
    my $self = shift;

    return unless $self->_verify_config_status();

    my @models = $self->models;
    return unless $self->_verify_models_are_assignable(@models);

    for my $model (@models) {
        $self->_assign_model($model);
    }

    $self->status_message('Models added to analysis project.');

    return 1;
}

sub _verify_config_status {
    my $self = shift;

    my $profile_item = $self->profile_item;
    if($profile_item->status eq 'active') {
        $self->error_message('Cannot manually assign a model to an active configuration.');
        return;
    }

    my $analysis_project = $profile_item->analysis_project;
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

    my $profile_item = $self->profile_item;
    my $analysis_project = $profile_item->analysis_project;

    my $config_text = Genome::Sys->read_file($profile_item->file_path);

    my $all_models_are_assignable = 1;
    for my $model (@models) {
        if($model->analysis_projects) {
            my $model_anp = $model->analysis_projects;
            unless($model_anp eq $analysis_project) {
                $self->error_message(
                    'Model %s is already assigned to another analysis project, %s.  Cannot assign it to %s.',
                    $model->id,
                    $model_anp->id,
                    $analysis_project->id,
                );
                $all_models_are_assignable = 0;
            }
        }

        my $model_class = $model->class;
        if(index($config_text, $model_class) < 0) {
            $self->error_message(
                'Model %s is not a type of model specified in the configuration profile item %s.  Cannot assign it to %s.',
                $model->id,
                $profile_item->id,
                $analysis_project->id,
            );
            $all_models_are_assignable = 0;
        }
    }

    return $all_models_are_assignable;
}

sub _assign_model {
    my $self = shift;
    my $model = shift;

    my $profile_item = $self->profile_item;
    my $analysis_project = $profile_item->analysis_project;

    if($model->analysis_projects) {
        $self->debug_message(
            'Model %s is already assigned to analysis project %s.',
            $model->id,
            $model->analysis_projects->id,
        );
        return 1;
    }

    $analysis_project->add_model_bridge(model => $model, config_profile_item => $profile_item);
    $self->debug_message(
        'Model %s assigned to analysis project %s.',
        $model->id,
        $model->analysis_projects->id,
    );

    return 1;
}

1;
