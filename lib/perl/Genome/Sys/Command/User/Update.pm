package Genome::Sys::Command::User::Update;

use strict;
use warnings;
use Genome;

class Genome::Sys::Command::User::Update {
    is => 'Command::V2',
    has => [
        user => {
            is => 'Genome::Sys::User',
            doc => 'user to modify',
        },
    ],
    has_optional => [
        username => {
            is => 'Text',
            doc => 'new username for user',
        },
        email => {
            is => 'Text',
            doc => 'new email for user',
        },
        name => {
            is => 'Text',
            doc => 'new name for user',
        },
    ],
    doc => 'update a user',
};

sub _is_hidden_in_docs {
    return !Genome::Sys->current_user_is_admin;
}

sub execute {
    my $self = shift;

    my $user = $self->user;

    if ($self->username) {
        $user->username($self->username);
    }

    if ($self->email) {
        $user->email($self->email);
    }

    if ($self->name) {
        $user->name($self->name);
    }

    $self->status_message("User successfully updated!");
    return 1;
}

1;

