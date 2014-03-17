package Genome::Model::Input;

use strict;
use warnings;

use Genome;

require Carp;

class Genome::Model::Input {
    table_name => 'model.model_input',
    type_name => 'genome model input',
    id_by => [
        value_class_name => {
            is => 'VARCHAR2',
            len => 255,
        },
        value_id => {
            is => 'VARCHAR2',
            len => 1000,
        },
        model_id => {
            is => 'Text',
            len => 32,
        },
        name => {
            is => 'VARCHAR2',
            len => 255,
        },
    ],
    has => [
        model => {
            is => 'Genome::Model',
            id_by => 'model_id',
            constraint_name => 'GMI_GM_FK',
        },
        value => {
            is => 'UR::Object',
            id_by => 'value_id',
            id_class_by => 'value_class_name',
        },
        filter_desc => {
            is => 'Text',
            valid_values => [ "forward-only", "reverse-only", undef ],
            is_optional => 1,
            doc => 'Filter to apply on the input value.',
        },
    ],
    has_optional => [
        _model_value => {
            is => 'Genome::Model',
            id_by => 'value_id',
        },
    ],
    schema_name => 'GMSchema',
    data_source => 'Genome::DataSource::GMSchema',
};

sub __display_name__ {
    my $self = shift;
    my $model = $self->model;
    my $value = $self->value;
    return (($model ? $model->__display_name__ : "") . " " . $self->name . ": " . ($value ? $value->__display_name__ : ""));
}

sub delete {
    my $self = shift;
    my $input_name = $self->__display_name__;

    # TODO I don't care for the instrument data special case. I think that inputs should be completely modifiable
    # without triggering builds being abandoned
    unless ($self->name eq 'instrument_data') {
        my $delete_rv = $self->SUPER::delete;
        Carp::confess "Could not delete input $input_name" unless $delete_rv;
        return 1;
    }

    for my $build ($self->builds_with_input) {
        $self->status_message("Abandoning build " . $build->__display_name__ . " that uses input $input_name");
        my $abandon_rv = eval { $build->abandon };
        unless ($abandon_rv) {
            Carp::confess "Could not abandon build " . $build->__display_name__ . " while deleting input";
        }
    }

    my $delete_rv = $self->SUPER::delete;
    Carp::confess "Could not delete input $input_name!" unless $delete_rv;
    return 1;
}

sub builds_with_input {
    my $self = shift;
    return unless $self->model;

    my @builds;
    for my $build ($self->model->builds) {
        next unless grep {
            $_->name eq $self->name and
            $_->value_id eq $self->value_id and
            $_->value_class_name eq $self->value_class_name
        } $build->inputs;
        push @builds, $build;
    }

    return @builds;
}

1;

