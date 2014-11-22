package Genome::ModelGroup::Command::Delete;

use strict;
use warnings;

use Genome;

class Genome::ModelGroup::Command::Delete {
    is => 'Command::V2',
    has => [
        model_group => {
            is => 'Genome::ModelGroup',
            shell_args_position => 1,
            require_user_verify => 1,
            doc => 'Model group to delete',
        },
    ]
};

sub help_synopsis {
    return <<"EOS"
    genome model-group delete --model-group 2
EOS
}

sub help_brief {
    return "delete a model-group";
}

sub help_detail {                           
    return <<EOS 
    delete a model-group
EOS
}

sub execute {
    my $self = shift;

    my $display_name = $self->model_group->__display_name__;
    unless ( $self->model_group->delete ) {
        $self->error_message("Cannot delete model-group: $display_name");
        return;
    }
    
    $self->status_message("Deleted model-group: $display_name");

    return 1;
}

1;
