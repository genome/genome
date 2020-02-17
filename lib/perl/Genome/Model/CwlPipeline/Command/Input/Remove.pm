package Genome::Model::CwlPipeline::Command::Input::Remove;

use strict;
use warnings;

use Genome;

class Genome::Model::CwlPipeline::Command::Input::Remove {
    is => 'Command::V2',
    has_input => [
        model => {
            is => 'Genome::Model::CwlPipeline',
            doc => 'the model(s) from which to remove inputs',
            is_many => 1,
            shell_args_position => 1,
        },
        name => {
            is => 'Text',
            doc => 'name of the input(s) to unset',
            shell_args_position => 2,
        },
        value => {
            is => 'Text',
            is_many => 1,
            doc => 'value(s) to set as input(s)',
            shell_args_position => 3,
        },
        value_class_name => {
            is => 'Text',
            doc => 'If set, the "value(s)" will be treated as IDs for objects of this class--resolution will be skipped',
            is_optional => 1,
        },
    ],
    doc => 'remove inputs from the model',
};

sub sub_command_category { 'input tools' }

sub help_detail {
    return <<EOHELP
Remove input values from a model.

For simple string inputs, the "value" should be the string itself.  For object inputs, the "value" should be the ID of the object to remove as an input.
EOHELP
}

sub execute {
    my $self = shift;

    my $count = 0;

    for my $model ($self->model) {
        VALUE: for my $value_identifier ($self->value) {

            my $value_class_name = $self->value_class_name;
            unless ($value_class_name) {
                my $value = $model->determine_input_object($self->name, $value_identifier);
                $value_class_name = $value->class;
                $value_identifier = $value->id;

                next VALUE unless $value;
            }

            my @input = Genome::Model::Input->get(
                model_id => $model->id,
                value_id => $value_identifier,
                value_class_name => $value_class_name,
                name => $self->name
            );

            for my $i (@input) {
                $i->delete;
                $count++;
            }
        }
    }

    $self->status_message('Removed %s input(s).', $count);
    return 1;
}

1;
