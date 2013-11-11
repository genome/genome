
package Genome::Model::Command::Remove;

use strict;
use warnings;

use Genome;
use Cwd;

class Genome::Model::Command::Remove {
    is => 'Genome::Command::Base',
    has => [
        models => {
            is => 'Genome::Model',
            shell_args_position => 1,
            is_many => 1,
            doc => 'The model to remove, specified by id, name or expression',
        },
        force_delete => {
            is => 'Boolean',
            default_value => 0,
            doc => 'A boolean flag to force model delete.(default_value=0)',
        },
    ],
    has_constant => [
        allow_model_with_builds => {
            is => 'Boolean',
            value => 0,
            doc => 'By default this command will not remove a model with builds, as this destroys records of past processing that has been done.',
        },
    ],
    doc => "delete a genome model, all of its builds, and logs",
};

sub help_synopsis {
    return <<"EOS"
genome model remove 12345
genome model remove mymodel 
genome model remove subject_name=FOO
EOS
}

sub _is_hidden_in_docs { return !Genome::Sys->current_user_is_admin }

sub execute {
    my $self = shift;

    my @models =  $self->models;
    my @names = map { $_->__display_name__ } @models;

    unless ($self->force_delete) {
        my $response = $self->_ask_user_question(
            'Are you sure you want to remove '
            . scalar(@models) 
            . ": @names?"
        );
        unless (defined $response and $response eq 'yes') {
            $self->status_message('Not deleting model(s).  Exiting.');
            return 1;
        }
    }

    for my $model (@models) {
        next if $model->isa("UR::DeletedRef");

        my @b = $model->builds;
        if(@b and not $self->allow_model_with_builds) {
            $self->error_message('Refusing to delete a model with builds');
            return;
        }

        $self->status_message("Removing model " . $model->__display_name__);
        my $model_id = $model->id;
        unless ($model->delete) {
            $self->error_message('Failed to delete model id '. $model_id);
            return;
        }
        $self->status_message('Succesfully removed model id '. $model_id);
    }

    return 1;
}




1;

