package Genome::Model::CwlPipeline::Command::Input::List;

use strict;
use warnings;

use Genome;

class Genome::Model::CwlPipeline::Command::Input::List {
    is => 'UR::Object::Command::List',
    has => [
        subject_class_name => {
            is_constant => 1,
            is => 'Text',
            value => 'Genome::Model::Input',
        },
    ],
    has_input => [
        model => {
            is => 'Genome::Model::CwlPipeline',
            doc => 'the model whose inputs to list',
        },
    ],
    has_param => [
        show => {
            is => 'Text',
            default => 'name,value',
        },
    ],
    doc => 'list the inputs for the model',
};

sub sub_command_category { 'input tools' }

sub _resolve_boolexpr {
    my $self = shift;

    return Genome::Model::Input->define_boolexpr(
        model_id => $self->model->id,
    );
}

1;
