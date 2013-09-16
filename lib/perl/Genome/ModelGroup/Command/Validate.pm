package Genome::ModelGroup::Command::Validate;

use strict;
use warnings;
use Genome;

class Genome::ModelGroup::Command::Validate {
    is => 'Genome::Command::Base',
    doc => "Validate a model group.",
    has_input => [
        model_group => {
            shell_args_position => 1,
            is => 'Genome::ModelGroup',
            doc => 'the model group to validate', 
        },
        distinct_subjects => {
            is => 'Boolean',
            is_optional => 1,
            doc => 'If true, make sure the models in this group have no duplicate subjects',
        },
    ],
};

sub execute {
    my $self = shift;

    my $failed = 0;
    if ($self->distinct_subjects) {
        $failed += $self->_has_duplicate_subjects;
    }

    if ($failed) {
        die $self->error_message("Model group " . $self->model_group->__display_name__ . " failed validation");
    } else {
        $self->status_message("Model group " . $self->model_group->__display_name__ . " succeeded validation");
    }

    return 1;
}

sub _has_duplicate_subjects {
    my $self = shift;

    my %models;
    for my $model ($self->model_group->models) {
        my $subject_id = $model->subject->id;
        push @{$models{$subject_id}}, $model->id;
    }

    my $failed = 0;
    for my $subject_id (keys %models) {
        my @models = @{$models{$subject_id}};
        if (scalar(@models) > 1) {
            $failed = 1;
            $self->error_message("Subject $subject_id occurs in more than one model: " . join(",", @models) );
        }
    }

    return $failed;
}

1;

