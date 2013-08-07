package Genome::Disk::Command::Allocation::Unarchive;

use strict;
use warnings;
use Genome;

class Genome::Disk::Command::Allocation::Unarchive {
    is => 'Genome::Disk::Command::Allocation::UnarchiveBase',
    has_optional => [
        allocations => {
            is => 'Genome::Disk::Allocation',
            is_many => 1,
            shell_args_position => 1,
            doc => 'allocations to be unarchived',
        },
        paths => {
            is => 'Text',
            is_many => 1,
            doc => 'pass a path instead of allocations and allocation will be looked up'
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

sub _execute {
    my $self = shift;

    if (!$self->allocations) {
        if ($self->paths) {
            my @allocs;
            for my $path ($self->paths) {
                my $alloc = Genome::Disk::Allocation->get_allocation_for_path($path);
                if ($alloc) { push @allocs, $alloc; }
            }
            $self->allocations(\@allocs);
        }

        if ($self->allocations) {
            $self->status_message('found allocation(s): ' . join("\n",map {$_->id} $self->allocations));
        } else {
            die 'You must provide either allocation id or the path you are trying to unarchive.';
        }
    }

    $self->status_message("Starting unarchive command...");

    for my $allocation ($self->allocations) {
        $self->debug_message("Unarchiving allocation " . $allocation->id);
        my $rv = $allocation->unarchive(reason => $self->reason);
        unless ($rv) {
            Carp::confess "Could not unarchive alloation " . $allocation->id;
        }
        $self->debug_message("Finished unarchiving allocation " . $allocation->id);
    }

    $self->status_message("Done unarchiving, exiting...");
    return 1;
}

1;
