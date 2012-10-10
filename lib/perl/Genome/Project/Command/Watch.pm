package Genome::Project::Command::Watch;

use strict;
use warnings;
use Genome;

class Genome::Project::Command::Watch {
    is => 'Command::V2',
    has => [
        project => {
            is => 'Genome::Project',
            doc => 'The project to watch.',
            shell_args_position => 1,
        },
        user => {
            is => 'Genome::Sys::User',
            is_optional => 1,
            doc => 'The user that will watch. If unspecified will be the current user.',
        },
    ],
};

sub execute {
    my $self = shift;

    unless ($self->user) {
        $self->user(Genome::Sys::User->get(username => Genome::Sys->username));
    }

    if (grep { $_ eq $self->user->id } $self->project->watcher_ids) {
        $self->status_message("User is already watching project.");
    } else {
        $self->status_message("Adding user as watcher.");
        $self->project->add_watcher($self->user);
    }

    return 1;
}

1;
