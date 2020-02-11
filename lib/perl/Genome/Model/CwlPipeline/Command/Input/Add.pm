package Genome::Model::CwlPipeline::Command::Input::Add;

use strict;
use warnings;

use Genome;

class Genome::Model::CwlPipeline::Command::Input::Add {
    is => 'Command::V2',
    has_input => [
        model => {
            is => 'Genome::Model::CwlPipeline',
            doc => 'the model(s) to which to add inputs',
            is_many => 1,
            shell_args_position => 1,
        },
        name => {
            is => 'Text',
            doc => 'name of the input(s) to set',
            shell_args_position => 2,
        },
        value => {
            is => 'Text',
            is_many => 1,
            doc => 'value(s) to set as input(s)',
            shell_args_position => 3,
        },
    ],
    doc => 'add inputs to the model',
};

sub sub_command_category { 'input tools' }

sub help_detail {
    return <<EOHELP
Add input values to a model.  Certain names will cause the value to be parsed into the relevant objects in the system; all other inputs will be added directly as strings.
EOHELP
}

sub execute {
    my $self = shift;

    for my $m ($self->model) {
        $m->process_input_data(
            $self->name, [$self->value]
        );
    }

    return 1;
}

1;
