package Genome::Sys::Command::User::Add;

use strict;
use warnings;
use Genome;
use Carp;

class Genome::Sys::Command::User::Add {
    is => 'Command::V2',
    has => [
        username => {
            is => 'Text',
            doc => 'system user name used for logging into systems (eg, jdoe)',
        },
        email => {
            is => 'Text',
            doc => 'email address at which the user can be contacted',
        },
    ],
    has_optional => [
        name => {
            is => 'Text',
            doc => 'name of the person this entry belongs to (eg, John Doe)',
        },
    ],
    doc => 'add a new user to the system',
};

sub _is_hidden_in_docs {
    return !Genome::Sys->current_user_is_admin;
}

sub execute {
    my $self = shift;

    my %params = (
        username => $self->username,
        email => $self->email,
    );
    $params{name} = $self->name if $self->name;
    
    my $user = Genome::Sys::User->create(%params);
    unless ($user) {
        Carp::confess "Could not create user!";
    }

    $self->status_message("Created user with username " . $self->username . " and email " . $self->email . "!");
    return 1;
}

1;

