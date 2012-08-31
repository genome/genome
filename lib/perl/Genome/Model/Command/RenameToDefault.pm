package Genome::Model::Command::RenameToDefault;

use strict;
use warnings;

use Genome;

class Genome::Model::Command::RenameToDefault {
    is => 'Genome::Command::Base',
    has => [
        models => {
            is => 'Genome::Model',
            is_many => 1,
            shell_args_position => 1,
            require_user_verify => 1,
            doc => 'The model(s) to rename to the default name.',
        },
    ],
    doc => "change the name of a model to the default",
};

sub help_synopsis {
    return;
}

sub execute {
    my $self = shift;

    for my $model ( $self->models ) {
        my $old_name = $model->name;
        $model->name('PLACEHOLDER'); # so the default name does not find this model
        $model->name( $model->default_model_name );
        $self->status_message("Rename $old_name to ".$model->name);
    }

    return 1;
}

1;

