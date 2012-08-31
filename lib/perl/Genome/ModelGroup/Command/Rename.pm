package Genome::ModelGroup::Command::Rename;

use strict;
use warnings;

use Genome;

class Genome::ModelGroup::Command::Rename {
    is => 'Genome::Command::Base',
    has => [
        from => {
            is => 'Genome::ModelGroup',
            id_by => 'group_id',
            shell_args_position => 1,
            doc => 'Current model group name or id.',
        },
        to => {
            is => 'Text',
            shell_args_position => 2,
            doc => 'New name for model group.',
        },
    ],
    doc => "change the name of a model-group (not the models in the group)",
};

sub sub_command_sort_position { 5 }

sub help_synopsis {
    return <<"EOS"
    genome model-group rename OLDNAME NEWNAME
    genome model-group rename 12345 NEWNAME
EOS
}

sub execute {
    my $self = shift;

    if ( $self->to eq $self->from->name ) {
        $self->error_message("New name is the same as the model's current name.");
        return;
    }

    my $old_name = $self->from->name;
    $self->from->rename( $self->to );

    unless ( $self->to eq $self->from->name ) {
        $self->error_message(
            sprintf(
                'Could not rename model group (<Id> %s <Name> %s) to new name (%s)', 
                $self->from->id, 
                $self->from->name, 
                $self->to,
            )
        );
        return;
    }

    return 1;
}

1;

