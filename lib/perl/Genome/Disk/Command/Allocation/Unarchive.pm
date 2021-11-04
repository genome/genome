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
        reason => {
            is => 'Text',
            is_optional => 1,
            doc => 'a note about why this unarchiving is taking place',
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

sub execute {
    my $self = shift;

    $self->status_message("Starting unarchive command...");

    for my $allocation ($self->allocations) {
        $self->status_message("Unarchiving allocation " . $allocation->id);
        my $rv = $allocation->unarchive(reason => $self->_reason);
        unless ($rv) {
            Carp::confess "Could not unarchive alloation " . $allocation->id;
        }
        $self->status_message("Finished unarchiving allocation " . $allocation->id);
    }

    $self->status_message("Done unarchiving, exiting...");
    return 1;
}

sub _reason {
    my $self = shift;

    my $reason = 'unarchived by ' . Genome::Sys->username;
    if (my $provided_reason = $self->reason) {
        $reason .= ' because: ' . $provided_reason;
    }

    return $reason;
}


1;
