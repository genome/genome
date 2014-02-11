package Genome::Project::Command::Update::Models;

use strict;
use warnings;

use Genome;
use Lingua::EN::Inflect;

class Genome::Project::Command::Update::Models {
    is => 'Command::V2',
    has_input => [
        'projects' => {
            is => 'Genome::Project',
            is_many => 1,
            shell_args_position => 1,
            doc => 'project(s) to be updated',
        },
        match_type => {
            is => 'Text',
            valid_values => ['instrument_data', 'sample'],
            default_value => 'instrument_data',
            doc => 'match models to the project based on these values',
        },
    ],
    has_constant => {
        model_user => {
            is => 'Text',
            doc => 'user whose models to add to the project',
            value => 'apipe-builder',
        },
    },
    doc => 'add models to a project that have matching instrument data or samples',
};

sub help_detail {
    return 'Add automatically-created models to a project that have instrument data or samples that are already parts of the project';
}

sub execute {
    my $self = shift;

    my $part_class = $self->_part_class_for_match_type($self->match_type);

    for my $project ($self->projects) {
        my @parts = $project->get_parts_of_class($part_class);

        my $model_getter = '_models_for_' . $self->match_type;
        my @models =
            grep { $_->created_by eq $self->model_user }
            $self->$model_getter(@parts);

        my @new_parts =
            map { $project->add_part(entity => $_) }
            grep { not $project->get_part($_) }
            @models;

        $self->status_message(
            'Added %s to project %s',
            Lingua::EN::Inflect::NO('model', scalar(@new_parts)),
            $project->__display_name__
        );
    }

    return 1;
}

my $MATCH_TYPE_TO_PART_CLASS = {
    instrument_data => 'Genome::InstrumentData',
    sample => 'Genome::Sample',
};

sub _part_class_for_match_type {
    my $self = shift;
    my $match_type = shift;

    return $MATCH_TYPE_TO_PART_CLASS->{$match_type};
}

sub _models_for_instrument_data {
    my $self = shift;
    my @parts = @_;

    my @inputs = Genome::Model::Input->get(value_id => [map $_->entity_id, @parts] );
    my @models = Genome::Model->get([map $_->model_id, @inputs]);
    return @models;
}

sub _models_for_sample {
    my $self = shift;
    my @parts = @_;

    my @models = Genome::Model->get(subject_id => [map $_->entity_id, @parts]);
    return @models;
}

1;
