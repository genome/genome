package Genome::Disk::Allocation;

use strict;
use warnings;

use Genome;
use Genome::Utility::Instrumentation;

use File::Copy::Recursive 'dircopy';
use Carp 'confess';

use List::Util 'shuffle';

our $TESTING_DISK_ALLOCATION = 0;

class Genome::Disk::Allocation {
    is => 'Genome::Notable',
    table_name => 'disk.allocation',
    id_by => [
        id => {
            is => 'Text',
            len => 64,
            doc => 'The id for the allocator event',
        },
    ],
    has => [
        disk_group_name => {
            is => 'Text',
            len => 40,
            doc => 'The name of the disk group',
        },
        mount_path => {
            is => 'Text',
            len => 255,
            doc => 'The mount path of the disk volume',
        },
        allocation_path => {
            is => 'Text',
            len => 4000,
            doc => 'The sub-dir of the disk volume for which space is allocated',
        },
        kilobytes_requested => {
            is => 'Number',
            len => 20,
            doc => 'The disk space allocated in kilobytes',
        },
        owner_class_name => {
            is => 'Text',
            len => 255,
            is_optional => 1,
            doc => 'The class name for the owner of this allocation',
        },
        owner_id => {
            is => 'Text',
            len => 255,
            is_optional => 1,
            doc => 'The id for the owner of this allocation',
        },
        owner => {
            is => 'UR::Object',
            id_by => 'owner_id',
            id_class_by => 'owner_class_name',
        },
        group_subdirectory => {
            is => 'Text',
            len => 255,
            doc => 'The group specific subdirectory where space is allocated',
        },
        absolute_path => {
            calculate_from => [ 'mount_path', 'group_subdirectory', 'allocation_path' ],
            calculate => q( $self->_absolute_path($mount_path, $group_subdirectory, $allocation_path); ),
        },
        volume => {
            is => 'Genome::Disk::Volume',
            calculate_from => 'mount_path',
            calculate => q( return Genome::Disk::Volume->get(mount_path => $mount_path, disk_status => 'active'); ),
        },
        group => {
            is => 'Genome::Disk::Group',
            calculate_from => 'disk_group_name',
            calculate => q( return Genome::Disk::Group->get(disk_group_name => $disk_group_name); ),
        },
        kilobytes_used_time => {
            is => 'DateTime',
            len => 11,
            is_optional => 1,
        },
        archive_after_time => {
            is => 'DateTime',
            len => 11,
            is_optional => 1,
            doc => 'After this time, this allocation is subject to being archived'
        },
    ],
    has_optional => [
        preserved => {
            is => 'Boolean',
            len => 1,
            default_value => 0,
            doc => 'If set, the allocation cannot be deallocated',
        },
        archivable => {
            is => 'Boolean',
            len => 1,
            default_value => 1,
            doc => 'If set, this allocation can be archived',
        },
        original_kilobytes_requested => {
            is => 'Number',
            len => 20,
            doc => 'The disk space allocated in kilobytes',
        },
        kilobytes_used => {
            is => 'Number',
            len => 20,
            default_value => 0,
            doc => 'The actual disk space used by owner',
        },
        creation_time => {
            is => 'DateTime',
            len => 11,
            doc => 'Time at which the allocation was created',
        },
        reallocation_time => {
            is => 'DateTime',
            len => 11,
            doc => 'The last time at which the allocation was reallocated',
        },
        owner_exists => {
            is => 'Boolean',
            calculate_from => [ 'owner_class_name', 'owner_id' ],
            calculate => q(
                my $owner_exists = eval { $owner_class_name->get($owner_id) };
                return $owner_exists ? 1 : 0;
            ),
        },
    ],
    schema_name => 'GMSchema',
    data_source => 'Genome::DataSource::GMSchema',
    id_generator => '-uuid',
};

# TODO This needs to be removed, site-specific
our @APIPE_DISK_GROUPS = qw/
    info_apipe
    info_apipe_ref
    info_alignments
    info_genome_models
    research
    systems_benchmarking
/;
our $CREATE_DUMMY_VOLUMES_FOR_TESTING = 1;
our @PATHS_TO_REMOVE; # Keeps track of paths created when no commit is on

sub _absolute_path {
    my $class = shift;
    my ($mount_path, $group_subdirectory, $allocation_path) = @_;
    return $mount_path .'/'. $group_subdirectory .'/'. $allocation_path;
}

sub du {
    my $self = shift;

    if ($self->is_archived) {
        die $self->error_message('Cannot `du` an archived allocation');
    }

    my $kb_used = 0;
    my $absolute_path = $self->absolute_path;
    if ( -d $absolute_path ) {
        # allow_errors will allow disk_usage_for_path to return a number even if
        # du emits errors (for instance, no read permissions for a subfolder)
        $kb_used = Genome::Sys->disk_usage_for_path(
            $absolute_path, allow_errors => 1) || 0;
    }

    return $kb_used;
}

sub allocate { return shift->create(@_); }
sub create {
    my ($class, %params) = @_;

    Genome::Utility::Instrumentation::inc('disk.allocation.create');

    # TODO Switch from %params to BoolExpr and pass in BX to autogenerate_new_object_id
    unless (exists $params{allocation_id}) {
        $params{allocation_id} = $class->__meta__->autogenerate_new_object_id;
    }

    # If no commit is on, make a dummy volume to allocate to
    if ($ENV{UR_DBI_NO_COMMIT}) {
        if ($CREATE_DUMMY_VOLUMES_FOR_TESTING) { # && !$params{mount_path}) {
            my $tmp_volume = Genome::Disk::Volume->create_dummy_volume(
                mount_path => $params{mount_path},
                disk_group_name => $params{disk_group_name},
            );
            $params{mount_path} = $tmp_volume->mount_path;
        }
    }

    my $self;
    Genome::Utility::Instrumentation::timer('disk.allocation.create', sub {
        $self = $class->_execute_system_command('_create', %params);
    });

    if ($self) {
        if ($ENV{UR_DBI_NO_COMMIT}) {
            push @PATHS_TO_REMOVE, $self->absolute_path;
        }
        else {
            $self->_log_change_for_rollback;
        }
    }


    return $self;
}

sub deallocate { return shift->delete(@_); }
sub delete {
    my ($class, %params) = @_;

    Genome::Utility::Instrumentation::inc('disk.allocation.delete');

    $class->_execute_system_command('_delete', %params);
    return 1;
}

sub reallocate {
    my ($class, %params) = @_;

    Genome::Utility::Instrumentation::inc('disk.allocation.reallocate');

    return $class->_execute_system_command('_reallocate', %params);
}

sub move {
    my ($class, %params) = @_;

    Genome::Utility::Instrumentation::inc('disk.allocation.move');

    return $class->_execute_system_command('_move', %params);
}

sub archive {
    my ($class, %params) = @_;

    Genome::Utility::Instrumentation::inc('disk.allocation.archive');

    return $class->_execute_system_command('_archive', %params);
}

sub unarchive {
    my ($class, %params) = @_;

    Genome::Utility::Instrumentation::inc('disk.allocation.unarchive');

    return $class->_execute_system_command('_unarchive', %params);
}

sub is_archived {
    my $self = shift;
    return $self->volume->is_archive;
}

sub tar_path {
    my $self = shift;
    return join('/', $self->absolute_path, 'archive.tar');
}

sub preserved {
    my ($self, $value, $reason) = @_;
    if (@_ > 1) {
        $reason ||= 'no reason given';
        $self->add_note(
            header_text => $value ? 'set to preserved' : 'set to unpreserved',
            body_text => $reason,
        );
        if ($value) {
            $self->_create_observer($self->_mark_read_only_closure($self->absolute_path));
        }
        else {
            $self->_create_observer($self->_set_default_permissions_closure($self->absolute_path));
        }
        return $self->__preserved($value);
    }
    return $self->__preserved();
}

sub archivable {
    my ($self, $value, $reason) = @_;
    if (@_ > 1) {
        $reason ||= 'no reason given';
        $self->add_note(
            header_text => $value ? 'set to archivable' : 'set to unarchivable',
            body_text => $reason,
        );
        return $self->__archivable($value);
    }
    return $self->__archivable;
}

sub _create {
    my $class = shift;
    my %params = @_;

    # Make sure that required parameters are provided
    my @missing_params;
    for my $param (qw/ disk_group_name allocation_path kilobytes_requested owner_class_name owner_id /) {
        unless (exists $params{$param} and defined $params{$param}) {
            push @missing_params, $param;
        }
    }
    if (@missing_params) {
        confess "Missing required params for allocation:\n" . join("\n", @missing_params);
    }

    # Make sure there aren't any extra params
    my $id = delete $params{allocation_id};
    $id = $class->__meta__->autogenerate_new_object_id unless defined $id; # TODO autogenerate_new_object_id should technically receive a BoolExpr
    my $kilobytes_requested = delete $params{kilobytes_requested};
    my $owner_class_name = delete $params{owner_class_name};
    my $owner_id = delete $params{owner_id};
    my $allocation_path = delete $params{allocation_path};
    my $disk_group_name = delete $params{disk_group_name};
    my $mount_path = delete $params{mount_path};
    my $exclude_mount_path = delete $params{exclude_mount_path};
    my $group_subdirectory = delete $params{group_subdirectory};
    my $kilobytes_used = delete $params{kilobytes_used} || 0;
    my $archive_after_time = delete $params{archive_after_time};
    
    if (%params) {
        confess "Extra parameters detected: " . Data::Dumper::Dumper(\%params);
    }

    unless ($owner_class_name->__meta__) {
        confess "Could not find meta information for owner class $owner_class_name, make sure this class exists!";
    }
    unless (defined $kilobytes_requested && $kilobytes_requested >= 0) {
        confess 'Kilobytes requested is not valid!';
    }
    if (my $parent_alloc = $class->get_parent_allocation($allocation_path)) {
        confess sprintf("Parent allocation (%s) found for %s", $parent_alloc->allocation_path, $allocation_path);
    }
    unless ($class->_verify_no_child_allocations($allocation_path)) {
        confess "Child allocation found for $allocation_path!";
    }
    unless (grep { $disk_group_name eq $_ } @APIPE_DISK_GROUPS) {
        confess "Can only allocate disk in apipe disk groups, not $disk_group_name. Apipe groups are: " . join(", ", @APIPE_DISK_GROUPS);
    }

    if ($ENV{GENOME_DB_PAUSE} and -e $ENV{GENOME_DB_PAUSE}) {
        print "Database updating has been paused; not going to attempt to allocate disk until the pause is released. Please stand by...\n";

        while (1) {
            sleep 30;
            last unless -e $ENV{GENOME_DB_PAUSE};
        }

        print "Database updating has been resumed, continuing allocation!\n";
    }

    my $group = Genome::Disk::Group->get(disk_group_name => $disk_group_name);
    confess "Could not find a group with name $disk_group_name" unless $group;
    if (defined $group_subdirectory and $group_subdirectory ne $group->subdirectory) {
        print STDERR "Given group subdirectory $group_subdirectory does not match retrieved group's subdirectory, ignoring provided value\n";
    }
    $group_subdirectory = $group->subdirectory;

    # If given a mount path, need to ensure it's valid by trying to get a disk volume with it. Also need to make
    # sure that the retrieved volume actually belongs to the supplied disk group and that it can be allocated to
    my @candidate_volumes;
    if (defined $mount_path) {
        $mount_path =~ s/\/$//; # mount paths in database don't have trailing /
        my $volume = Genome::Disk::Volume->get(mount_path => $mount_path, disk_status => 'active', can_allocate => 1);
        confess "Could not get volume with mount path $mount_path" unless $volume;

        unless (grep { $_ eq $disk_group_name } $volume->disk_group_names) {
            confess "Volume with mount path $mount_path is not in supplied group $disk_group_name!";
        }

        my @reasons;
        push @reasons, 'disk is not active' if $volume->disk_status ne 'active';
        push @reasons, 'allocation turned off for this disk' if $volume->can_allocate != 1;

        if ($exclude_mount_path && $volume->mount_path eq $exclude_mount_path) {
            push @reasons, 'Specified mount path matched the excluded mount path.';
        }

        if (@reasons) {
            confess "Requested volume with mount path $mount_path cannot be allocated to:\n" . join("\n", @reasons);
        }

        push @candidate_volumes, $volume;
    } else {
        my %candidate_volume_params = (
            disk_group_name => $disk_group_name,
        );
        if (defined $exclude_mount_path) {
            $candidate_volume_params{'exclude'} = $exclude_mount_path;
        }
        push @candidate_volumes, $class->_get_candidate_volumes(
            %candidate_volume_params);
    }

    my %parameters = (
        disk_group_name              => $disk_group_name,
        kilobytes_requested          => $kilobytes_requested,
        original_kilobytes_requested => $kilobytes_requested,
        allocation_path              => $allocation_path,
        owner_class_name             => $owner_class_name,
        owner_id                     => $owner_id,
        group_subdirectory           => $group_subdirectory,
        id                           => $id,
        creation_time                => UR::Context->current->now,
        archive_after_time           => $archive_after_time,        
    );

    # Make sure that we never attempt to create an allocation that has an absolute path that already exists. There are
    # several places in this module that may delete a newly-created "candidate" allocation under the assumption that the
    # absolute path for that allocation is empty. If there are preexisting files they will be unintentionally deleted
    # when the candidate allocation is destroyed.
    #
    # For example: a user is trying to create an allocation for a specific path that already exists without an
    # allocation (by specifying mount_path, disk_group_name, and allocation_path, they can force a specific absolute
    # path for the new allocation). If a candidate allocation is created at this path and then destroyed for some
    # arbitrary reason, the user will lose their files.

    for my $candidate_volume (@candidate_volumes) {
        my $candidate_path = $class->_absolute_path($candidate_volume->mount_path, $group_subdirectory, $allocation_path);
        if ( -e $candidate_path ) {
            confess "The allocation path $candidate_path already exists. If you are attempting to create an allocation "
                . "for an existing path, please move the path to a temporary location before continuing.";
        }
    }

    my $self = $class->_get_allocation_without_lock(\@candidate_volumes, \%parameters);

    $self->debug_message(sprintf("Allocation (%s) created at %s",
        $id, $self->absolute_path));

    # a restrictive umask can break builds for other users; force it to be friendly
    umask(0002);
    # If we cannot create the directory delete the new allocation
    my $dir = eval {
        Genome::Sys->create_directory($self->absolute_path);
    };
    my $error = $@;
    unless (defined($dir) and ( -d $dir ) and not $error) {
        $class->error_message(sprintf(
                "Failed to create directory (%s) with return value = '%s', and error:\n%s",
                $self->absolute_path, $dir || '', $error));
        $self->delete;
        confess $error;
    }

    return $self;
}

sub _get_allocation_without_lock {
    my ($class, $candidate_volumes, $parameters) = @_;
    my $kilobytes_requested = $parameters->{'kilobytes_requested'};

    # We randomize to avoid the rare repeated contention case
    my @randomized_candidate_volumes = (@$candidate_volumes,
        shuffle(@$candidate_volumes));

    my $chosen_allocation;
    for my $candidate_volume (@randomized_candidate_volumes) {
        if ($candidate_volume->allocated_kb + $kilobytes_requested
                <= $candidate_volume->soft_limit_kb) {
            my $candidate_allocation = $class->SUPER::create(
                mount_path => $candidate_volume->mount_path,
                %$parameters,
            );
            unless ($candidate_allocation) {
                die 'Failed to create candidate allocation';
            }
            _commit_unless_testing();

            # Reload so we guarantee that we calculate the correct allocated_kb
            if (not $ENV{UR_DBI_NO_COMMIT}) {
                UR::Context->current->reload($candidate_volume);
            }

            if ($candidate_volume->is_allocated_over_soft_limit) {
                $class->status_message(sprintf("%s's allocated_kb exceeded soft limit (%d kB), rolling back allocation.", $candidate_volume->mount_path, $candidate_volume->soft_limit_kb, 'kB'));
                $candidate_allocation->delete();
                _commit_unless_testing();
            } elsif ($candidate_volume->is_used_over_soft_limit) {
                $class->status_message(sprintf("%s's used_kb exceeded soft limit (%d %s), rolling back allocation.", $candidate_volume->mount_path, $candidate_volume->soft_limit_kb, 'kB'));
                $candidate_allocation->delete();
                _commit_unless_testing();
            } else {
                $chosen_allocation = $candidate_allocation;
                last;
            }
        }
    }

    unless (defined $chosen_allocation) {
        Carp::confess $class->error_message(sprintf(
            "Could not create allocation in specified disk group (%s), which contains %d volumes:\n%s\n",
            $parameters->{disk_group_name}, scalar(@$candidate_volumes),
            join("\n", map { $_->mount_path } @$candidate_volumes),
        ));
    }

    return $chosen_allocation;
}

sub _delete {
    my ($class, %params) = @_;
    my $id = delete $params{allocation_id};
    confess "Require allocation ID!" unless defined $id;
    if (%params) {
        confess "Extra params found: " . Data::Dumper::Dumper(\%params);
    }

    my $self = $class->get($id);
    if ($self->preserved) {
        confess sprintf("Allocation (%s) has been marked as preserved, cannot deallocate!",
            $self->id);
    }

    my $path = $self->absolute_path;

    $self->SUPER::delete;

    $class->_create_observer(
        $class->_mark_for_deletion_closure($path),
        $class->_remove_directory_closure($path),
    );

    return 1;
}

sub _reallocate {
    my ($class, %params) = @_;
    my $id = delete $params{allocation_id};
    confess "Require allocation ID!" unless defined $id;
    my $kilobytes_requested = delete $params{kilobytes_requested} || 0;
    my $allow_reallocate_with_move = delete $params{allow_reallocate_with_move};
    my $grow_only = delete $params{grow_only};
    if (%params) {
        confess "Found extra params: " . Data::Dumper::Dumper(\%params);
    }

    my $mode = $class->_retrieve_mode();
    my $self = $class->$mode($id);
    my $old_kb_requested = $self->kilobytes_requested;

    my $kb_used = $self->du();

    # Cache kilobytes used
    $self->kilobytes_used($kb_used);

    my $actual_kb_requested = List::Util::max($kb_used, $kilobytes_requested);
    if ($grow_only && ($actual_kb_requested <= $old_kb_requested)) {
        $self->status_message(
            "Not changing kilobytes_requested, because grow_only = 1 & actual usage < original_kilobytes_requested");
        return 1;
    }
    if ($actual_kb_requested > $kilobytes_requested) {
        $self->status_message(sprintf(
                "Setting kilobytes_requested to %s based on `du` for allocation %s",
                $actual_kb_requested, $self->id));
    }

    $self->reallocation_time(UR::Context->current->now);
    $self->kilobytes_requested($actual_kb_requested);
    _commit_unless_testing();

    my $volume = $self->volume;
    my $succeeded;
    if ($volume->allocated_kb < $volume->hard_limit_kb) {
        $succeeded = 1;
    } else {
        if ($allow_reallocate_with_move) {
            my $old_mount_path = $self->mount_path;
            $self->move();
            unless ($self->mount_path eq $old_mount_path) {
                $succeeded = 1;
            } else {
                $self->error_message(sprintf(
                        "Failed to reallocate and move allocation %s from volume %s.",
                        $self->id, $volume->mount_path));
            }
        } else {
            $self->error_message(sprintf(
                    "Failed to reallocate allocation %s on volume %s because the volume is beyond its quota.",
                    $self->id, $volume->mount_path));
        }
    }

    unless ($succeeded) {
        # Rollback kilobytes_requested
        my $max_kilobytes_requested = List::Util::max($kb_used, $old_kb_requested);
        my $msg = $old_kb_requested == $max_kilobytes_requested ? 'Rolling back' : 'Setting';

        $self->status_message(sprintf(
                "%s kilobytes_requested to %d for allocation %s.",
                $msg, $max_kilobytes_requested, $self->id));

        $self->kilobytes_requested($max_kilobytes_requested);
        _commit_unless_testing();
    }

    return $succeeded;
}

sub move_shadow_path {
    my $allocation_path = shift;
    return sprintf("%s-move_allocation_destination", $allocation_path)
}

sub _move {
    my ($class, %params) = @_;
    my $id = delete $params{allocation_id};

    my $kilobytes_requested = delete $params{kilobytes_requested};
    my $group_name = delete $params{disk_group_name};
    my $new_mount_path = delete $params{target_mount_path};

    if (%params) {
        confess "Extra parameters given to allocation move method: " . join(',', sort keys %params);
    }

    # get allocation lock
    my $allocation_lock = Genome::Disk::Allocation->get_lock($id);
    unless ($allocation_lock) {
        confess "Could not lock allocation with ID $id";
    }
    my $mode = $class->_retrieve_mode;
    my $self = Genome::Disk::Allocation->$mode($id);
    unless ($self) {
        Genome::Sys->unlock_resource(resource_lock => $allocation_lock);
        confess "Found no allocation with ID $id";
    }

    my $original_absolute_path = $self->absolute_path;

    # make shadow allocation
    my %creation_params = (
        disk_group_name => $self->disk_group_name,
        kilobytes_requested => $self->kilobytes_requested,
        owner_class_name => "UR::Value",
        owner_id => "shadow_allocation",
        exclude_mount_path => $self->mount_path,
        allocation_path => move_shadow_path($self->allocation_path),
    );

    # I think that it's dangerous to specify the new mount path, but this
    # feature existed, so nnutter and I kept it during this refactor.
    if ($new_mount_path) {
        $creation_params{'mount_path'} = $new_mount_path;
    }

    if ($group_name) {
        $creation_params{disk_group_name} = $group_name;
    }

    # The shadow allocation is just a way of keeping track of our temporary
    # additional disk usage during the move.
    my $shadow_allocation = $class->create(%creation_params);
    my $shadow_absolute_path = $shadow_allocation->absolute_path;

    # copy files to shadow allocation
    my $copy_rv = eval {
        Genome::Sys->rsync_directory(
            source_directory => $original_absolute_path,
            target_directory => $shadow_absolute_path,
        );
    };
    unless ($copy_rv) {
        my $copy_error_message = $@;
        if (-d $shadow_absolute_path) {
            if (Genome::Sys->remove_directory_tree($shadow_absolute_path)) {
                $shadow_allocation->delete;
            }
        }
        Genome::Sys->unlock_resource(resource_lock => $allocation_lock);
        confess(sprintf(
                "Could not copy allocation %s from %s to %s: %s",
                $self->id, $original_absolute_path, $shadow_absolute_path, $copy_error_message));
    }

    my $new_volume_final_path = $class->_absolute_path(
        $shadow_allocation->mount_path,
        $shadow_allocation->group_subdirectory,
        $self->allocation_path);


    Genome::Sys->create_directory($new_volume_final_path);
    unless (rename $shadow_allocation->absolute_path, $new_volume_final_path) {
        Genome::Sys->unlock_resource(resource_lock => $allocation_lock);
        my $shadow_allocation_abs_path = $shadow_allocation->absolute_path;
        $shadow_allocation->delete;
        confess($self->error_message(sprintf(
                "Could not move shadow allocation path (%s) to final path (%s).  This should never happen, even when 100%% full.",
                $shadow_allocation_abs_path, $self->absolute_path)));
    }

    # Change the shadow allocation to reserve some disk on the old volume until
    # those files are deleted.
    my $old_mount_path = $self->mount_path;
    $self->mount_path($shadow_allocation->mount_path);
    $shadow_allocation->mount_path($old_mount_path);

    # This is here because certain objects (Build & SoftwareResult) don't
    # calculate their data_directories from their disk_allocations.
    $self->_update_owner_for_move;

    _symlink($self->absolute_path, $original_absolute_path);

    _commit_unless_testing();

    Genome::Sys->unlock_resource(resource_lock => $allocation_lock);

    $shadow_allocation->delete;

    $class->_create_observer(
        $class->_mark_for_deletion_closure($original_absolute_path),
        $class->_remove_directory_closure($original_absolute_path),
    );
}

sub _archive {
    my ($class, %params) = @_;
    my $id = delete $params{allocation_id};
    if (%params) {
        confess "Extra parameters given to allocation move method: " . join(',', sort keys %params);
    }

    # Lock and load allocation object
    my $allocation_lock = Genome::Disk::Allocation->get_lock($id);
    unless ($allocation_lock) {
        confess "Could not lock allocation with ID $id";
    }
    my $mode = $class->_retrieve_mode;
    my $self = Genome::Disk::Allocation->$mode($id);
    unless ($self) {
        Genome::Sys->unlock_resource(resource_lock => $allocation_lock);
        confess "Found no allocation with ID $id";
    }

    my $current_allocation_path = $self->absolute_path;
    my $archive_allocation_path = join('/', $self->volume->archive_mount_path, $self->group_subdirectory, $self->allocation_path);
    my $tar_path = join('/', $archive_allocation_path, 'archive.tar');

    # This gets set to true immediately before tarball creation is started. This allows for conditional clean up of the
    # archive directory in case of failure, which is nice because that requires an LSF job be scheduled.
    my $tarball_created = 0;

    eval {
        if ($self->is_archived) {
            confess "Allocation " . $self->id . " is already archived!";
        }
        if (!$self->archivable) {
            confess "Allocation " . $self->id . " is not flagged as archivable!";
        }
        unless (-e $current_allocation_path) {
            confess "Allocation path $current_allocation_path does not exist!";
        }

        unless (Genome::Sys->recursively_validate_directory_for_read_write_access($current_allocation_path)) {
            confess "Some files in allocation directory $current_allocation_path are not readable or writable, this must be fixed before archiving can continue";
        }

        # Reallocate so we're reflecting the correct size at time of archive.
        $self->reallocate();
        if (!$ENV{UR_DBI_NO_COMMIT} and $self->kilobytes_requested < 1048576) { # Must be greater than 1GB to be archived
            confess(sprintf(
                "Total size of files at path %s is only %s, which is not greater than 1GB!",
                $current_allocation_path, $self->kilobytes_requested));
        }

        my $mkdir_cmd = "mkdir -p $archive_allocation_path";
        my $cd_cmd = "cd $current_allocation_path";
        my $tar_cmd = "/bin/ls -A | tar --create --file $tar_path -T -";
        my $cmd = join(' && ', $mkdir_cmd, $cd_cmd, $tar_cmd);

        $tarball_created = 1;
        # It's very possible that if no commit is on, the volumes/allocations being dealt with are test objects that don't
        # exist out of this local UR context, so bsubbing jobs would fail.
        if ($ENV{UR_DBI_NO_COMMIT}) {
            Genome::Sys->shellcmd(cmd => $cmd);
        }
        else {
            my ($job_id, $status);

            my @signals = qw/ INT TERM /;
            for my $signal (@signals) {
                $SIG{$signal} = sub {
                    print STDERR "Cleanup activated within allocation, cleaning up LSF jobs\n";
                    eval { Genome::Sys->kill_lsf_job($job_id) } if $job_id;
                    $self->_cleanup_archive_directory($archive_allocation_path);
                    die "Received signal, exiting.";
                }
            }

            ($job_id, $status) = Genome::Sys->bsub_and_wait(
                queue => $ENV{GENOME_ARCHIVE_LSF_QUEUE},
                job_group => '/archive',
                log_file => "\"/tmp/$id\"", # Entire path must be wrapped in quotes because older allocation IDs contain spaces
                cmd => "\"$cmd\"",          # If the command isn't wrapped in quotes, the '&&' is misinterpreted by
                                            # bash (rather than being "bsub '1 && 2' it is looked at as 'bsub 1' && '2')
            );

            for my $signal (@signals) {
                delete $SIG{$signal};
            }

            unless ($status eq 'DONE') {
                confess "LSF job $job_id failed to execute $cmd, exited with status $status";
            }
        }

        $self->mount_path($self->volume->archive_mount_path);
        $self->_update_owner_for_move;

        my $rv = $self->_commit_unless_testing;
        confess "Could not commit!" unless $rv;
    };
    my $error = $@; # Record error so it can be investigated after unlocking

    # If only there were finally blocks...
    Genome::Sys->unlock_resource(resource_lock => $allocation_lock) if $allocation_lock;

    if ($error) {
        eval {
            $self->_cleanup_archive_directory($archive_allocation_path) if $tarball_created;
        };
        my $cleanup_error = $@;
        my $msg = "Could not archive allocation " . $self->id . ", error:\n$error";
        if ($cleanup_error) {
            $msg .= "\n\nWhile cleaning up archive tarball, encoutered error:\n$cleanup_error";
        }
        confess $msg;
    }

    # Never make filesystem changes if no commit is enabled
    unless ($ENV{UR_DBI_NO_COMMIT}) {
        my $rv = Genome::Sys->remove_directory_tree($current_allocation_path);
        unless ($rv) {
            confess "Could not remove unarchived allocation path $current_allocation_path. " .
                "Database changes have been committed, so clean this up manually!";
        }
    }

    return 1;
}

sub unarchive_shadow_path {
    my $allocation_path = shift;
    return sprintf("%s-unarchive_allocation_destination", $allocation_path);
}

sub _unarchive {
    my ($class, %params) = @_;
    my $id = delete $params{allocation_id};
    my $reason = delete $params{reason};
    if (%params) {
        confess "Extra parameters given to allocation unarchive method: " . join(',', sort keys %params);
    }

    # Lock and load allocation object
    my $allocation_lock = Genome::Disk::Allocation->get_lock($id);
    unless ($allocation_lock) {
        confess "Could not lock allocation with ID $id";
    }
    my $mode = $class->_retrieve_mode;
    my $self = Genome::Disk::Allocation->$mode($id);
    unless ($self) {
        Genome::Sys->unlock_resource(resource_lock => $allocation_lock);
        confess "Found no allocation with ID $id";
    }

    unless ($self->is_archived) {
        Genome::Sys->unlock_resource(resource_lock => $allocation_lock);
        $self->status_message("Allocation is not archived, cannot unarchive. Exiting.");
        return 1;
    }

    my %creation_params = (
        disk_group_name => $self->disk_group_name,
        kilobytes_requested => $self->kilobytes_requested,
        owner_class_name => "UR::Value",
        owner_id => "shadow_allocation",
        allocation_path => unarchive_shadow_path($self->allocation_path),
    );
    # shadow_allocation ensures that we wont over allocate our destination volume
    my $shadow_allocation = $class->create(%creation_params);
    unless ($shadow_allocation) {
        confess(sprintf("Failed to create shadow allocation from params %s",
                Data::Dumper::Dumper(\%creation_params)));
    }

    my $archive_path = $self->absolute_path;
    my $target_path = $shadow_allocation->absolute_path;

    # Inferred path prior to archiving so we can symlink the new allocation,
    # to this old location.
    (my $old_absolute_path = $archive_path) =~ s/^\/gscarchive/\/gscmnt/;

    my $tar_path = $self->tar_path;
    my $cmd = "tar -C $target_path -xf $tar_path";

    eval {
        # It's very possible that if no commit is on, the volumes/allocations being dealt with are test objects that don't
        # exist out of this local UR context, so bsubbing jobs would fail.
        if ($ENV{UR_DBI_NO_COMMIT}) {
            Genome::Sys->shellcmd(cmd => $cmd);
        } else {
            # If this process should be killed, the LSF job needs to be cleaned up
            my ($job_id, $status);

            # Signal handlers are added like this so an anonymous sub can be used, which handles variables defined
            # in an outer scope differently than named subs (in this case, $job_id, $self, and $archive_path).
            my @signals = qw/ INT TERM /;
            for my $signal (@signals) {
                $SIG{$signal} = sub {
                    print STDERR "Cleanup activated within allocation, cleaning up LSF jobs\n";
                    eval { Genome::Sys->kill_lsf_job($job_id) } if $job_id;
                    $self->_cleanup_archive_directory($archive_path);
                    die "Received signal, exiting.";
                }
            }

            ($job_id, $status) = Genome::Sys->bsub_and_wait(
                queue => $ENV{GENOME_ARCHIVE_LSF_QUEUE},
                job_group => '/unarchive',
                log_file => "\"/tmp/$id\"", # Entire path must be wrapped in quotes because older allocation IDs contain spaces
                cmd => "\"$cmd\"", # If the command isn't wrapped in quotes, the '&&' is misinterpreted by
                                   # bash (rather than being "bsub '1 && 2' it is looked at as 'bsub 1' && '2')
            );

            for my $signal (@signals) {
                delete $SIG{$signal};
            }

            unless ($status eq 'DONE') {
                confess "Could not execute command $cmd via LSF job $job_id, received status $status";
            }
        }

        # Make updates to the allocation
        $self->mount_path($shadow_allocation->volume->mount_path);
        Genome::Sys->create_directory($self->absolute_path);
        unless (rename $shadow_allocation->absolute_path, $self->absolute_path) {
            confess($self->error_message(sprintf(
                    "Could not move shadow allocation path (%s) to final path (%s).  This should never happen, even when 100%% full.",
                    $shadow_allocation->absolute_path, $self->absolute_path)));
        }
        $self->_update_owner_for_move;
        $self->archivable(0, 'allocation was unarchived'); # Wouldn't want this to be immediately re-archived... trolololol

        if ($old_absolute_path ne $self->absolute_path) {
            _symlink($old_absolute_path, $self->absolute_path);
        }

        unless ($self->_commit_unless_testing) {
            confess "Could not commit!";
        }
    };
    my $error = $@;

    # finally blocks would be really sweet. Alas...
    Genome::Sys->unlock_resource(resource_lock => $allocation_lock) if $allocation_lock;
    $shadow_allocation->delete();

    if ($error) {
        if ($target_path and -d $target_path and not $ENV{UR_DBI_NO_COMMIT}) {
            Genome::Sys->remove_directory_tree($target_path);
        }
        confess "Could not unarchive, received error:\n$error";
    } else {
        $self->add_note(
            header_text => 'unarchived',
            body_text => $reason || 'no reason given',
        );
    }

    $self->_cleanup_archive_directory($archive_path);
    return 1;
}

# Locks the allocation, if lock is not manually released (it had better be!) it'll be automatically
# cleaned up on program exit
sub get_lock {
    my ($class, $id, $tries) = @_;
    $tries ||= 60;
    my $allocation_lock = Genome::Sys->lock_resource(
        resource_lock => $ENV{GENOME_LOCK_DIR} . '/allocation/allocation_' . join('_', split(' ', $id)),
        max_try => $tries,
        block_sleep => 1,
    );
    # lock implies mutation so object should be reloaded
    $class->_reload_allocation($id);
    return $allocation_lock;
}

sub has_valid_owner {
    my $self = shift;
    my $meta = $self->owner_class_name->__meta__;
    return 0 unless $meta;
    return 1;
}

sub __display_name__ {
    my $self = shift;
    return $self->absolute_path;
}

# Using a system call when not in dev mode is a hack to get around the fact that we don't
# have software transactions. Allocation need to be able to make its changes and commit
# immediately so locks can be released in a timely manner. Without software transactions,
# the only two ways of doing this are to: just commit away, or create a subprocess and commit
# there. Comitting in the calling process makes it possible for objects in an intermediate
# state to be committed unintentionally (for example, if a user makes an object of type Foo,
# then creates an allocation for it, and then intends to finish instantiating it after),
# which can either lead to outright rejection by the database if a constraint is violated or
# other more subtle problems in the software. So, that leaves making a subprocess, which is
# slow but won't lead to other problems.
our @_execute_system_command_perl5opt = '-MGenome';
sub _execute_system_command {
    my ($class, $method, %params) = @_;
    if (ref($class)) {
        $params{allocation_id} = $class->id;
        $class = ref($class);
    }
    confess "Require allocation ID!" unless exists $params{allocation_id};

    my $allocation;
    if ($ENV{UR_DBI_NO_COMMIT}) {
        $allocation = $class->$method(%params);
    }
    else {
        # remove the parens if you DARE
        my @includes = map { ( '-I' => $_ ) } UR::Util::used_libs;

        my $param_string = Genome::Utility::Text::hash_to_string(\%params, 'q');
        my $perl_program_string = sprintf("%s->%s(%s); UR::Context->commit;",
            $class, $method, $param_string);

        my @cmd = (
            'genome-perl',
            @includes,
            @_execute_system_command_perl5opt,
            '-e',
            $perl_program_string
        );
        unless (system(@cmd) == 0) {
            my $msg = "Could not perform allocation action!";
            if ($@) {
                $msg .= " Error: $@";
            }
            confess $msg;
        }
        $allocation = $class->_reload_allocation($params{allocation_id});
    }


    return $allocation;
}

sub _log_change_for_rollback {
    my $self = shift;
    # If the owner gets rolled back, then delete the allocation. Make sure the allocation hasn't already been deleted,
    # which can happen if the owner is coded well and cleans up its own mess during rollback.
    my $remove_sub = sub {
        my $allocation_id = $self->id;
        $self->unload;
        my $loaded_allocation = Genome::Disk::Allocation->get($allocation_id);
        $loaded_allocation->delete if ($loaded_allocation);
    };
    my $allocation_change = UR::Context::Transaction->log_change(
        $self->owner, 'UR::Value', $self->id, 'external_change', $remove_sub,
    );
    return 1;
}

# Some owners track their absolute path separately from the allocation, which means they also need to be
# updated when the allocation is moved. That special logic goes here
sub _update_owner_for_move {
    my $self = shift;
    my $owner_class = $self->owner_class_name;
    eval "require $owner_class";
    my $error = $@;
    if ($error) {
        if ($error =~ /Can't locate [\w|\/|\.]+ in \@INC/) {
            # Some allocations are owned by classes that no longer exist. In this case,
            # we just return successfully, as there's no owner to update.
            return 1;
        }
        else {
            die "Could not load owner class $owner_class, received error: $error";
        }
    }

    my $owner = $self->owner;
    return 1 unless $owner;

    if ($owner->isa('Genome::SoftwareResult')) {
        $owner->output_dir($self->absolute_path);
    }
    elsif ($owner->isa('Genome::Model::Build')) {
        $owner->data_directory($self->absolute_path);
    }

    return 1;
}

# Unloads the allocation and then reloads to ensure that changes from database are retrieved
sub _reload_allocation {
    my ($class, $id) = @_;

    my $allocation;
    my $mode = $class->_retrieve_mode;
    if ($mode eq 'get') {
        $allocation = Genome::Disk::Allocation->get($id);
    } elsif ($mode eq 'load') {
        $allocation = UR::Context->current->reload($class, id => $id);
    } else {
        die 'Unrecognized _retrieve_mode: ' . $class->_retrieve_mode;
    }

    return $allocation;
}

# Creates an observer that executes the supplied closures
sub _create_observer {
    my ($class, @closures) = @_;
    my $observer;
    my $callback = sub {
        $observer->delete if $observer;
        for my $closure (@closures) {
            &$closure;
        }
    };

    if ($ENV{UR_DBI_NO_COMMIT}) {
        &$callback;
        return 1;
    }

    $observer = UR::Context->add_observer(
        aspect => 'commit',
        callback => $callback,
    );
    return 1;
}

# Returns a closure that removes the given locks
sub _unlock_closure {
    my ($class, @locks) = @_;
    return sub {
        for my $lock (@locks) {
            Genome::Sys->unlock_resource(resource_lock => $lock) if -e $lock;
        }
    };
}

# Returns a closure that creates a directory at the given path
sub _create_directory_closure {
    my ($class, $path) = @_;
    return sub {
        # This method currently returns the path if it already exists instead of failing
        my $dir = eval{ Genome::Sys->create_directory($path) };
        if (defined $dir and -d $dir) {
            chmod(02775, $dir);
        }
        else {
            print STDERR "Could not create allocation directcory at $path!\n";
            print "$@\n" if $@;
        }
    };
}

# Returns a closure that removes the given directory
sub _remove_directory_closure {
    my ($class, $path) = @_;
    return sub {
        if (-d $path and not $ENV{UR_DBI_NO_COMMIT}) {
            print STDERR "Removing allocation directory $path\n";
            my $rv = Genome::Sys->remove_directory_tree($path);
            unless (defined $rv and $rv == 1) {
                Carp::cluck "Could not remove allocation directory $path!";
            }
        }
    };
}

# Make a file at the root of the allocation directory indicating that the allocation is gone,
# which makes it possible to figure out which directories should have been deleted but failed.
sub _mark_for_deletion_closure {
    my ($class, $path) = @_;
    return sub {
        # In some instances -d $path will return true even if $path was
        # recently moved (or possibly deleted?). This causes a "touch
        # $path/ALLOCATION_DELETED" even though $path does not exist, causing
        # touch to fail.  Attempting to opendir $path (even though it doesn't
        # exist) before testing -d $path seems to fix the issue. This is
        # probably an issue with NFS caching everything in 'struct stat', which
        # is used by stat(), which is used by -d. We tried to stat($path) and
        # system("stat $path") before testing -d, and in both cases stat()
        # returned information about a non-existent path.
        # See here: http://irccrew.org/~cras/nfs-coding-howto.html#attrcache (visited Nov 27 2012)
        my $rv = opendir(my $dh, $path);
        closedir $dh if ($rv);
        if (-d $path and not $ENV{UR_DBI_NO_COMMIT}) {
            print STDERR "Marking directory at $path as deallocated\n";
            system("touch $path/ALLOCATION_DELETED");
        }
    };
}

# Changes an allocation directory to read-only
sub _mark_read_only_closure {
    my ($class, $path) = @_;
    return sub {
        return unless -d $path and not $ENV{UR_DBI_NO_COMMIT};

        require File::Find;
        sub mark_read_only {
            my $file = $File::Find::name;
            if (-d $file) {
                chmod 0555, $file;
            }
            else {
                chmod 0444, $file
            }
        };

        print STDERR "Marking directory at $path read-only\n";
        File::Find::find(\&mark_read_only, $path);
    };
}

# Changes an allocation directory to default permissions
sub _set_default_permissions_closure {
    my ($class, $path) = @_;
    return sub {
        return unless -d $path and not $ENV{UR_DBI_NO_COMMIT};

        require File::Find;
        sub set_default_perms {
            my $file = $File::Find::name;
            if (-d $file) {
                chmod 0775, $file;
            }
            else {
                chmod 0664, $file
            }
        };

        print STDERR "Setting permissions to defaults for $path\n";
        File::Find::find(\&set_default_perms, $path);
    };
}

# Class method for determining if the given path has a parent allocation
sub _verify_no_parent_allocation {
    my ($class, $path) = @_;
    unless (defined $path) {
        Carp::croak("no path for parent check");
    }
    my $allocation = $class->get_parent_allocation($path);
    return !(defined $allocation);
}

# Returns parent allocation for the given path if one exists
sub get_parent_allocation {
    my ($class, $path) = @_;
    Carp::confess("no path defined") unless defined $path;
    my ($allocation) = $class->get(allocation_path => $path);
    return $allocation if $allocation;

    my $dir = File::Basename::dirname($path);
    if ($dir ne '.' and $dir ne '/') {
        return $class->get_parent_allocation($dir);
    }
    return;
}

sub _allocation_path_from_full_path {
    my ($class, $path) = @_;
    my $allocation_path = $path;
    my $mount_path = $class->_get_mount_path_from_full_path($path);
    return unless $mount_path;

    my $group_subdir = $class->_get_group_subdir_from_full_path_and_mount_path($path, $mount_path);
    return unless $group_subdir;

    $allocation_path =~ s/^$mount_path//;
    $allocation_path =~ s/^\/$group_subdir//;
    $allocation_path =~ s/^\///;
    return $allocation_path;
}

sub _get_mount_path_from_full_path {
    my ($class, $path) = @_;
    my @parts = grep { defined $_ and $_ ne '' } split(/\//, $path);
    for (my $i = 0; $i < @parts; $i++) {
        my $volume_subpath = '/' . join('/', @parts[0..$i]);
        my ($volume) = Genome::Disk::Volume->get(mount_path => $volume_subpath);
        return $volume_subpath if $volume;
    }
    return;
}

sub _get_group_subdir_from_full_path_and_mount_path {
    my ($class, $path, $mount_path) = @_;
    my $subpath = $path;
    $subpath =~ s/$mount_path//;
    $subpath =~ s/\///;
    my @parts = split(/\//, $subpath);

    for (my $i = 0; $i < @parts; $i++) {
        my $group_subpath = join('/', @parts[0..$i]);
        my ($group) = Genome::Disk::Group->get(subdirectory => $group_subpath);
        return $group_subpath if $group;
    }
    return;
}

# Checks for allocations beneath this one, which is also invalid
sub _verify_no_child_allocations {
    my ($class, $path) = @_;

    $path =~ s/\/+$//;

    my $meta        = $class->__meta__;
    my $table_name  = $meta->table_name;
    my $data_source = $meta->data_source;
    my $owner       = $data_source->owner;
    my $query_string;

    if ($data_source->isa('UR::DataSource::Oracle')) {
        my $fq_table_name = join('.', $owner, $table_name);
        $query_string = sprintf(q(select * from %s where allocation_path like ? AND rownum <= 1), $fq_table_name);
    } elsif ($data_source->isa('UR::DataSource::Pg') || $data_source->isa('UR::DataSource::SQLite')) {
        $query_string = sprintf(q(select * from %s where allocation_path like ? LIMIT 1), $table_name);
    } else {
        $class->error_message("Falling back on old child allocation detection behavior.");
        return !($class->et_child_allocations($path));
    }

    my $dbh = $data_source->get_default_handle();
    my $query_object = $dbh->prepare($query_string);
    $query_object->bind_param(1, $path . "/%");
    $query_object->execute();

    my $row_arrayref = $query_object->fetchrow_arrayref();
    my $err = $dbh->err;
    my $errstr = $dbh->errstr;
    $query_object->finish();

    if ($err) {
        die $class->error_message(sprintf(
                "Could not verify no child allocations: %s", $errstr));
    }

    return !defined $row_arrayref;
}

sub get_all_allocations_for_path {
    my ($class, $path) = @_;
    return ($class->get_parent_allocation($path), $class->get_child_allocations($path))
}

sub get_child_allocations {
    my ($class, $path) = @_;
    $path =~ s/\/+$//;
    return $class->get('allocation_path like' => $path . '/%');
}

# Makes sure the supplied kb amount is valid (nonzero and bigger than mininum)
sub _check_kb_requested {
    my ($class, $kb) = @_;
    return 0 unless defined $kb;
    return 1;
}

# Returns a list of volumes that meets the given criteria
sub _get_candidate_volumes {
    my ($class, %params) = @_;

    my $disk_group_name = delete $params{disk_group_name};
    my $exclude = delete $params{exclude};

    if (%params) {
        confess "Illegal arguments to _get_candidate_volumes: " . join(', ', keys %params);
    }

    my %volume_params = (
        disk_group_names => $disk_group_name,
        can_allocate => 1,
        disk_status => 'active',
    );

    # 'not like' caused conversion error on Oracle but 'not in' with anonymous array works
    $volume_params{'mount_path not in'} = [$exclude] if $exclude;
    # XXX Shouldn't need this, 'archive' should obviously be a status.
    #       This might be a performance issue.
    my @volumes = grep { not $_->is_archive } Genome::Disk::Volume->get(
        %volume_params, '-order_by' => ['-cached_unallocated_kb']);
    unless (@volumes) {
        confess "Did not get any allocatable and active volumes belonging to group $disk_group_name.";
    }

    return @volumes;
}


# When no commit is on, ordinarily an allocation goes to a dummy volume that only exists locally. Trying to load
# that dummy volume would lead to an error, so use a get instead.
sub _retrieve_mode {
    return 'get' if $ENV{UR_DBI_NO_COMMIT};
    return 'load';
}

sub _cleanup_archive_directory {
    my ($class, $directory) = @_;
    my $cmd = "if [ -a $directory] ; then rm -rf $directory ; fi";
    unless ($ENV{UR_DBI_NO_COMMIT}) {
        my ($job_id, $status) = Genome::Sys->bsub_and_wait(
            queue => $ENV{GENOME_ARCHIVE_LSF_QUEUE},
            cmd => "\"$cmd\"",
        );
        confess "Failed to execute $cmd via LSF job $job_id, received status $status" unless $status eq 'DONE';
    }
    return 1;
}

# Cleans up directories, useful when no commit is on and the test doesn't clean up its allocation directories
# or in the case of reallocate with move when a copy fails and temp data needs to be removed
sub remove_test_paths {
    for my $path (@PATHS_TO_REMOVE) {
        next unless -d $path;
        Genome::Sys->remove_directory_tree($path);
        if ($ENV{UR_DBI_NO_COMMIT}) {
            __PACKAGE__->debug_message("Removing allocation path $path because UR_DBI_NO_COMMIT is on");
        }
        else {
            __PACKAGE__->debug_message("Cleaning up allocation path $path");
        }
    }
}
END {
    if (@PATHS_TO_REMOVE) {
        remove_test_paths();
    }
}

sub get_allocation_for_path {
    my ($class, $path) = @_;

    my @parts = split(/\//, $path);
    @parts = @parts[4..$#parts]; # Remove mount path and group subdirectory

    my $allocation;
    # Try finding allocation by allocation path, removing subdirectories from the end after each attempt
    while (@parts) {
        $allocation = Genome::Disk::Allocation->get(allocation_path => join('/', @parts));
        last if $allocation;
        @parts = @parts[0..($#parts - 1)];
    }

    return $allocation;
}

sub _symlink {
    # to fix previously existing/broken symlinks by restoring chain
    my $real_path = shift;
    my $old_path = shift;

    my $old_parent_dir = File::Basename::dirname($old_path);
    if(! -e $old_parent_dir){
        Genome::Sys->create_directory($old_parent_dir);
    }
    if(-d $old_parent_dir && ! -e $old_path) {
        symlink $real_path, $old_path;
    }

    return (-l $old_path && readlink($old_path) eq $real_path);
}

sub _commit_unless_testing {
    if ($TESTING_DISK_ALLOCATION || !$ENV{UR_DBI_NO_COMMIT}) {
        return UR::Context->commit();
    } else {
        return 1;
    }
}

1;
