
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
    doc => "delete a genome model, all of its builds, and logs",
};

sub help_synopsis {
    return <<"EOS"
genome model remove 12345
genome model remove mymodel 
genome model remove subject_name=FOO
EOS
}

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

