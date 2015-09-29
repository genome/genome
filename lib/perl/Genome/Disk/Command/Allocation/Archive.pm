package Genome::Disk::Command::Allocation::Archive;

use strict;
use warnings;
use Genome;
use Carp;

class Genome::Disk::Command::Allocation::Archive {
    is => 'Command::V2',
    has_optional => [
        allocations => {
            is => 'Genome::Disk::Allocation',
            is_many => 1,
            shell_args_position => 1,
            doc => 'allocations to be archived',
        },
        paths => {
            is => 'Text',
            is_many => 1,
            doc => 'paths for which to look up the corresponding allocations to be archived',
        },
    ],
    doc => 'archives the given allocations',
};

sub help_detail {
    return 'archives the given allocations, putting them on tape and removing them from filesystem';
}

sub help_brief {
    return 'archives the given allocations';
}

sub execute {
    my $self = shift;
    $self->status_message("Starting archive command...");

    my @allocations = $self->allocations;
    if ($self->paths) {
        for my $path ($self->paths) {
            my $allocation = Genome::Disk::Allocation->get_allocation_for_path($path);
            if ($allocation) {
                $self->debug_message('Found allocation %s for path: %s', $allocation->id, $path);
                push @allocations, $allocation;
            } else {
                die $self->error_message('No allocation found for path: %s', $path);
            }
        }
    }

    for my $allocation (@allocations) {
        $self->status_message("Archiving allocation " . $allocation->id);
        if ($allocation->archivable == 0) {
            $self->status_message("Skipping allocation " . $allocation->id . ", not set to archivable");
            next;
        }

        my $rv = $allocation->archive;
        unless ($rv) {
            Carp::confess "Could not archive allocation " . $allocation->id;
        }
        $self->status_message("Finished archiving allocation " . $allocation->id);
    }

    $self->status_message("Done archiving, exiting...");
    return 1;
}

1;

