package Genome::Disk::Command::Allocation::Purge;

use strict;
use warnings;
use Genome;
use Carp;

class Genome::Disk::Command::Allocation::Purge {
    is => 'Command::V2',
    has => [
        allocations => {
            is => 'Genome::Disk::Allocation',
            is_many => 1,
            shell_args_position => 1,
            doc => 'allocations to be purged',
        },
        reason => {
            is => 'Text',
            doc => 'why these allocations are being purged',
        },
    ],
    doc => 'purges the given allocations',
};

sub help_detail {
    return "PERMANENTLY DELETES the allocation's files. The files are "
        . 'temporarily stored in a trash folder for short-term recovery. '
        . 'Database records of this allocation will be kept.';
}

sub help_brief {
    return 'PERMANENTLY DELETES files';
}

sub execute {
    my $self = shift;
    $self->status_message("Starting purge command...");

    for my $allocation ($self->allocations) {
        if ($allocation->status eq 'purged') {
            $self->status_message("Skipping allocation %s, "
                . "status already 'purged'", $allocation->id);
            next;
        }
        $self->status_message("Purging allocation %s...", $allocation->id);

        my $rv = $allocation->purge(reason => $self->reason);
        unless ($rv) {
            Carp::confess "Could not purge allocation " . $allocation->id;
        }
        $self->status_message("Finished purging allocation %s.",
            $allocation->id);
    }

    $self->status_message("Done purging, exiting...");
    return 1;
}

1;
