package Genome::Disk::Command::Allocation::Unarchive;

use strict;
use warnings;
use Genome;

class Genome::Disk::Command::Allocation::Unarchive {
    is => 'Command::V2',
    has => [
        allocations => {
            is => 'Genome::Disk::Allocation',
            is_many => 1,
            shell_args_position => 1,
            doc => 'allocations to be unarchived',
        },
    ],
    doc => 'unarchives the given allocations',
};

sub help_detail {
    return 'unarchives the given allocations, moving them from tape back onto the filesystem';
}

sub help_brief {
    return 'unarchives the given allocations';
}

sub _is_hidden_in_docs {
    return !Genome::Sys->current_user_has_role('archive');
}

sub execute {
    my $self = shift;
    $self->status_message("Starting unarchive command...");

    for my $allocation ($self->allocations) {
        $self->debug_message("Unarchiving allocation " . $allocation->id);
        my $rv = $allocation->unarchive;
        unless ($rv) {
            Carp::confess "Could not unarchive alloation " . $allocation->id;
        }
        $self->debug_message("Finished unarchiving allocation " . $allocation->id);
    }

    $self->status_message("Done unarchiving, exiting...");
    return 1;
}

1;

