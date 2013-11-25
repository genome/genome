package Genome::Disk::Detail::Allocation::Creator;

use strict;
use warnings;

use Genome;

use Carp qw(confess);
use List::Util 'shuffle';


class Genome::Disk::Detail::Allocation::Creator {
    has => {
        parameters => {
            is => 'Genome::Disk::Detail::Allocation::CreationParameters',
        },
    },
};


sub create_allocation {
    my $self = shift;

    $self->verify_no_parent_or_child_allocations;
    $self->wait_for_database_pause;

    my @candidate_volumes = $self->candidate_volumes;

    my $allocation_object = $self->_get_allocation_without_lock(
        \@candidate_volumes);

    $self->create_directory_or_delete_allocation($allocation_object);

    Genome::Timeline::Event::Allocation->created('initial creation',
        $allocation_object);

    return $allocation_object;
}

sub verify_no_parent_or_child_allocations {
    my $self = shift;

    $self->verify_no_parent_allocation;
    $self->verify_no_child_allocations;
}

sub verify_no_parent_allocation {
    my $self = shift;

    if (my $parent_alloc = Genome::Disk::Allocation->get_parent_allocation(
            $self->parameters->allocation_path)) {
        confess sprintf("Parent allocation (%s) found for %s",
            $parent_alloc->allocation_path, $self->parameters->allocation_path);
    }
}

sub verify_no_child_allocations {
    my $self = shift;

    unless (Genome::Disk::Allocation->_verify_no_child_allocations(
            $self->parameters->allocation_path)) {
        confess sprintf("Child allocation found for %s!",
            $self->parameters->allocation_path);
    }
}

sub wait_for_database_pause {
    if ($ENV{GENOME_DB_PAUSE} and -e $ENV{GENOME_DB_PAUSE}) {
        print "Database updating has been paused; not going to attempt "
            . "to allocate disk until the pause is released. "
            . "Please stand by...\n";

        while (1) {
            sleep 30;
            last unless -e $ENV{GENOME_DB_PAUSE};
        }

        print "Database updating has been resumed, continuing allocation!\n";
    }
}

sub candidate_volumes {
    my $self = shift;

    my @candidate_volumes;
    Genome::Utility::Instrumentation::timer(
        'disk.allocation.create.candidate_volumes.selection', sub {
            if (defined $self->parameters->mount_path) {
                @candidate_volumes = $self->candidate_volumes_from_mount_path;
            } else {
                @candidate_volumes =
                    $self->candidate_volumes_without_mount_path;
            }
    });

    return @candidate_volumes;
}

sub candidate_volumes_from_mount_path {
    my $self = shift;

    my $mount_path = $self->parameters->mount_path;
    my $disk_group_name = $self->parameters->disk_group_name;

    my $volume = Genome::Disk::Volume->get(mount_path => $mount_path,
        disk_status => 'active', can_allocate => 1);
    confess "Could not get volume with mount path $mount_path" unless $volume;

    unless (grep { $_ eq $disk_group_name } $volume->disk_group_names) {
        confess sprintf(
            "Volume with mount path %s is not in supplied group %s!",
            $mount_path, $disk_group_name);
    }

    return ($volume);
}

sub candidate_volumes_without_mount_path {
    my $self = shift;

    my %candidate_volume_params = (
        disk_group_name => $self->parameters->disk_group_name,
    );

    my $exclude_mount_path = $self->parameters->exclude_mount_path;
    if (defined $exclude_mount_path) {
        $candidate_volume_params{'exclude'} = $exclude_mount_path;
    }
    return $self->_get_candidate_volumes(
        %candidate_volume_params);
}

# Returns a list of volumes that meets the given criteria
sub _get_candidate_volumes {
    my ($self, %params) = @_;

    my $disk_group_name = delete $params{disk_group_name};
    my $exclude = delete $params{exclude};

    if (%params) {
        confess "Illegal arguments to _get_candidate_volumes: "
            . join(', ', keys %params);
    }

    my %volume_params = (
        disk_group_names => $disk_group_name,
        can_allocate => 1,
        disk_status => 'active',
    );

    # 'not like' caused conversion error on Oracle,
    # but 'not in' with anonymous array works
    $volume_params{'mount_path not in'} = [$exclude] if $exclude;
    # XXX Shouldn't need this, 'archive' should obviously be a status.
    #       This might be a performance issue.
    my @volumes = grep { not $_->is_archive } Genome::Disk::Volume->get(
        %volume_params, '-order_by' => ['-cached_unallocated_kb']);
    unless (@volumes) {
        confess sprintf(
            "Did not get any allocatable, active volumes from group %s.",
            $disk_group_name);
    }

    return @volumes;
}

sub _get_allocation_without_lock {
    my ($self, $candidate_volumes) = @_;

    my $allocation_object;
    Genome::Utility::Instrumentation::timer(
        'disk.allocation.create.get_allocation_without_lock', sub {
            $allocation_object = $self->_get_allocation_without_lock_impl(
                $candidate_volumes);
    });
    return $allocation_object;
}

sub _get_allocation_without_lock_impl {
    my ($self, $candidate_volumes) = @_;
    # We randomize to avoid the rare repeated contention case
    my @randomized_candidate_volumes = (@$candidate_volumes,
        shuffle(@$candidate_volumes));

    my $chosen_allocation;
    for my $candidate_volume (@randomized_candidate_volumes) {
        if ($candidate_volume->has_space(
                $self->parameters->kilobytes_requested)) {
            $self->_verify_allocation_path_unused($candidate_volume);

            $chosen_allocation = $self->_attempt_allocation_creation(
                $candidate_volume);
            if ($chosen_allocation) {
                last;
            }
        }
    }

    unless (defined $chosen_allocation) {
        Carp::confess $self->error_message(sprintf(
            "Could not create allocation in specified disk group (%s), "
            . "which contains %d volumes:\n%s\n",
            $self->parameters->disk_group_name, scalar(@$candidate_volumes),
            join("\n", map { $_->mount_path } @$candidate_volumes),
        ));
    }

    $chosen_allocation->debug_message(sprintf("Allocation (%s) created at %s",
        $chosen_allocation->id, $chosen_allocation->absolute_path));

    return $chosen_allocation;
}

sub _verify_allocation_path_unused {
    my ($self, $candidate_volume) = @_;

    # Make sure that we never attempt to create an allocation that has an
    # absolute path that already exists. There are several places in this
    # module that may delete a newly-created "candidate" allocation under the
    # assumption that the absolute path for that allocation is empty. If there
    # are preexisting files they will be unintentionally deleted when the
    # candidate allocation is destroyed.
    #
    # For example: a user is trying to create an allocation for a specific path
    # that already exists without an allocation (by specifying mount_path,
    # disk_group_name, and allocation_path, they can force a specific absolute
    # path for the new allocation). If a candidate allocation is created at
    # this path and then destroyed for some arbitrary reason, the user will
    # lose their files.

    Genome::Utility::Instrumentation::timer(
        'disk.allocation.create.candidate_volumes'
        .  '.existing_allocation_path_check', sub {
            my $candidate_path = Genome::Disk::Allocation->_absolute_path(
                $candidate_volume->mount_path,
                $self->parameters->group_subdirectory,
                $self->parameters->allocation_path);
            if ( -e $candidate_path ) {
                confess sprintf("The allocation path %s already exists. "
                    . "If you are attempting to create an allocation "
                    . "for an existing path, please move the path to a "
                    .  "temporary location before continuing.",
                    $candidate_path);
            }
    });
}

sub _attempt_allocation_creation {
    my ($self, $candidate_volume) = @_;

    my $candidate_allocation = Genome::Disk::Allocation->SUPER::create(
        mount_path => $candidate_volume->mount_path,
        $self->parameters->as_hash,
    );
    unless ($candidate_allocation) {
        die 'Failed to create candidate allocation';
    }
    Genome::Disk::Allocation::_commit_unless_testing();

    # Reload so we guarantee that we calculate the correct allocated_kb
    if (not $ENV{UR_DBI_NO_COMMIT}) {
        UR::Context->current->reload($candidate_volume);
    }

    if ($candidate_volume->is_allocated_over_soft_limit) {
        Genome::Utility::Instrumentation::inc('disk.allocation."
            . "get_allocation_without_lock.rollback.over_allocated');
        $self->status_message(sprintf(
                "%s's allocated_kb exceeded soft limit (%d kB), "
                . "rolling back allocation.",
                $candidate_volume->mount_path,
                $candidate_volume->soft_limit_kb, 'kB'));
        $candidate_allocation->delete();
        Genome::Disk::Allocation::_commit_unless_testing();
        return;

    } elsif ($candidate_volume->is_used_over_soft_limit) {
        Genome::Utility::Instrumentation::inc('disk.allocation."
            . "get_allocation_without_lock.rollback.over_used');
        $self->status_message(sprintf(
                "%s's used_kb exceeded soft limit (%d %s), "
                . "rolling back allocation.",
                $candidate_volume->mount_path,
                $candidate_volume->soft_limit_kb, 'kB'));
        $candidate_allocation->delete();
        Genome::Disk::Allocation::_commit_unless_testing();
        return;
    }

    return $candidate_allocation;
}

sub create_directory_or_delete_allocation {
    my ($self, $allocation_object) = @_;

    # If we cannot create the directory delete the new allocation
    my $dir;
    eval {
        Genome::Utility::Instrumentation::timer(
            'disk.allocation.create.create_directory', sub {
                $dir = Genome::Sys->create_directory(
                    $allocation_object->absolute_path);
        });
    };
    my $error = $@;
    unless (defined($dir) and ( -d $dir ) and not $error) {
        $self->error_message(sprintf(
                "Failed to create directory (%s) with return value = '%s', "
                . "and error:\n%s", $allocation_object->absolute_path,
                $dir || '', $error));
        $allocation_object->delete;
        confess $error;
    }
}


1;
