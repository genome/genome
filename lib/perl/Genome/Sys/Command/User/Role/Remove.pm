package Genome::Sys::Command::User::Role::Remove;

use strict;
use warnings;
use Genome;

class Genome::Sys::Command::User::Role::Remove {
    is => 'Command::V2',
    has => [
        roles => {
            is => 'Genome::Sys::User::Role',
            shell_args_position => 1,
            is_many => 1,
            doc => 'roles to be removed',
        },
    ],
    doc => 'removes role(s) from the system',
};

sub _is_hidden_in_docs {
    return !Genome::Sys->current_user_is_admin;
}

sub execute {
    my $self = shift;

    my @names = map { $_->name } $self->roles;
    
    for my $role ($self->roles) {
        my $name = $role->name;
        my $rv = $role->delete;
        unless ($rv) {
            Carp::confess "Could not remove role $name!";
        }
    }

    $self->status_message("Successfuly removed roles: " . join(', ', @names));
    return 1;
}


1;

