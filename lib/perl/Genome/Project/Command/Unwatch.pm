package Genome::Project::Command::Unwatch;

use strict;
use warnings;
use Genome;

class Genome::Project::Command::Unwatch {
    is => 'Command::V2',
    has => [
        project => {
            is => 'Genome::Project',
            doc => 'The project to unwatch.',
            shell_args_position => 1,
        },
        user => {
            is => 'Genome::Sys::User',
            is_optional => 1,
            doc => 'The user that will unwatch. If unspecified will be the current user.',
        },
    ],
};

sub execute {
    my $self = shift;

    unless ($self->user) {
        $self->user(Genome::Sys::User->get(username => Genome::Sys->username));
    }

    if (grep { $_ eq $self->user->id } $self->project->watcher_ids) {
        $self->status_message("Removing user as watcher.");
        $self->project->remove_watcher($self->user);
    } else {
        $self->status_message("User is not watching project.");
    }

    return 1;
}

1;
