package Genome::Disk::Allocation;

use strict;
use warnings;

use Genome;
use Genome::Utility::Instrumentation;

use Carp qw(croak confess);
use Digest::MD5 qw(md5_hex);
use File::Copy::Recursive qw(dircopy dirmove);
use File::Find;
use File::Find::Rule;
use Cwd;
use DateTime;

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
        status => {
            is => 'Text',
            valid_values => ['active', 'completed', 'purged', 'archived', 'invalid'],
            doc => 'The current status of this allocation',
            default_value => 'active',
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
            default_value => &_default_archive_after_time,
            is_optional => 1,
            doc => 'After this time, this allocation is subject to being archived'
        },
        timeline_events => {
            is_many => 1,
            is => 'Genome::Timeline::Event::Allocation',
            reverse_as => 'object',
        },
    ],
    has_optional => [
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
        file_summaries => {
            is => 'Genome::Disk::Allocation::FileSummary',
            is_many => 1,
            reverse_as => 'allocation'
        }
    ],
    schema_name => 'GMSchema',
    data_source => 'Genome::DataSource::GMSchema',
    id_generator => '-uuid',
};

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


sub set_permissions_read_only {
    my $self = shift;

    $self->set_file_permissions(0444);
    $self->set_directory_permissions(0555);

    chmod 0555, $self->absolute_path;
}

sub set_file_permissions {
    my ($self, $mode) = @_;
    my @files = File::Find::Rule->file->in($self->absolute_path);
    chmod $mode, @files;
}

sub set_directory_permissions {
    my ($self, $mode) = @_;
    my @subdirs = File::Find::Rule->directory->in($self->absolute_path);
    chmod $mode, @subdirs;
}


sub allocate { return shift->create(@_); }
sub create {
    my ($class, %params) = @_;

    Genome::Utility::Instrumentation::inc('disk.allocation.create');

    # TODO Switch from %params to BoolExpr and pass in BX to autogenerate_new_object_id
    unless (exists $params{id}) {
        $params{id} = $class->__meta__->autogenerate_new_object_id;
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

sub archivable {
    my ($self, $value, $reason) = @_;
    if (@_ > 1) {
        $reason ||= 'no reason given';
        $self->add_note(
            header_text => $value ? 'set to archivable' : 'set to unarchivable',
            body_text => $reason,
        );
        my $event = $value ? 'unpreserved' : 'preserved';
        Genome::Timeline::Event::Allocation->$event($reason, $self);
        return $self->__archivable($value);
    }
    return $self->__archivable;
}

sub _create {
    my $class = shift;

    my $pars = Genome::Disk::Detail::Allocation::CreationParameters->create(@_);

    my $creator = Genome::Disk::Detail::Allocation::Creator->create(
        parameters => $pars);
    return $creator->create_allocation();
}


sub _delete {
    my ($class, %params) = @_;
    my $id = delete $params{id};
    confess "Require allocation ID!" unless defined $id;
    if (%params) {
        confess "Extra params found: " . Data::Dumper::Dumper(\%params);
    }

    my $self = $class->get($id);

    my $path = $self->absolute_path;

    for my $event ($self->timeline_events) {
        $event->delete();
    }

    $self->SUPER::delete;

    $class->_create_observer(
        $class->_mark_for_deletion_closure($path),
        $class->_remove_directory_closure($path),
    );

    return 1;
}

sub _reallocate {
    my $class = shift;

    my %parameters = @_;
    $parameters{allocation_id} = delete $parameters{id};

    my $reallocator = Genome::Disk::Detail::Allocation::Reallocator->create(
        %parameters);

    return $reallocator->reallocate;
}

sub move_shadow_path {
    my $allocation_path = shift;
    return sprintf("%s-move_allocation_destination", $allocation_path)
}

sub move_shadow_params {
    my $self = shift;
    return (
        disk_group_name => $self->disk_group_name,
        kilobytes_requested => $self->kilobytes_requested,
        owner_class_name => "UR::Value",
        owner_id => "shadow_allocation",
        exclude_mount_path => $self->mount_path,
        allocation_path => move_shadow_path($self->allocation_path),
    );
}

sub _move {
    my $class = shift;

    my %parameters = @_;
    $parameters{allocation_id} = delete $parameters{id};

    my $mover = Genome::Disk::Detail::Allocation::Mover->create(
        %parameters);

    return $mover->move;
}

sub _archive {
    my ($class, %params) = @_;
    my $id = delete $params{id};
    if (%params) {
        confess "Extra parameters given to allocation move method: " . join(',', sort keys %params);
    }

    my ($self, $allocation_lock) = $class->get_with_lock($id);

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
    } else {
        Genome::Timeline::Event::Allocation->archived(
            'archived',
            $self,
        )
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

sub unarchive_shadow_params {
    my $self = shift;
    return (
        disk_group_name => $self->disk_group_name,
        kilobytes_requested => $self->kilobytes_requested,
        owner_class_name => "UR::Value",
        owner_id => "shadow_allocation",
        allocation_path => unarchive_shadow_path($self->allocation_path),
    );
}

sub _unarchive {
    my ($class, %params) = @_;
    my $id = delete $params{id};
    my $reason = (delete $params{reason} || 'no reason given');
    if (%params) {
        confess "Extra parameters given to allocation unarchive method: " . join(',', sort keys %params);
    }

    my ($self, $allocation_lock) = $class->get_with_lock($id);

    unless ($self->is_archived) {
        Genome::Sys->unlock_resource(resource_lock => $allocation_lock);
        $self->status_message("Allocation is not archived, cannot unarchive. Exiting.");
        return 1;
    }

    my %creation_params = $self->unarchive_shadow_params;
    # shadow_allocation ensures that we wont over allocate our destination volume
    my $shadow_allocation = $class->shadow_get_or_create(%creation_params);

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
        $self->archive_after_time(Genome::Disk::Command::Allocation::DelayArchiving->_resolve_date_from_months(3));

        if ($old_absolute_path ne $self->absolute_path) {
            _symlink_new_path_from_old($old_absolute_path, $self->absolute_path);
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
            body_text => $reason,
        );
        Genome::Timeline::Event::Allocation->unarchived($reason, $self);
    }

    $self->_cleanup_archive_directory($archive_path);
    return 1;
}

# Locks the allocation, if lock is not manually released (it had better be!)
# it'll be automatically cleaned up on program exit
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

sub get_with_lock {
    my $class = shift;
    my $id = shift;

    my $lock = $class->get_lock($id);
    unless ($lock) {
        confess "Could not lock allocation with ID $id";
    }
    my $mode = $class->_retrieve_mode;
    my $self = Genome::Disk::Allocation->$mode($id);
    unless ($self) {
        Genome::Sys->unlock_resource(resource_lock => $lock);
        confess "Found no allocation with ID $id";
    }

    return ($self, $lock);
}

sub shadow_get_or_create {
    my $class = shift;
    my %params = @_;

    my %create_extra_params;
    for my $param (qw(exclude_mount_path kilobytes_requested)) {
        if ($params{$param}) {
            $create_extra_params{$param} = delete $params{$param};
        }
    }

    unless ($params{allocation_path}) {
        croak 'allocation_path is required for shadow_get_or_create';
    }

    my $md5_hex = md5_hex($params{allocation_path});;

    my $path = File::Spec->join(
        $ENV{GENOME_LOCK_DIR},
        'allocation',
        'allocation_path_' . $md5_hex,
    );

    my $lock = Genome::Sys->lock_resource(
        resource_lock => $path,
        wait_announce_interval => 600,
    );

    my $allocation = $class->get(%params);
    if ($allocation) {
        if ($create_extra_params{kilobytes_requested}) {
            $allocation->reallocate(
                kilobytes_requested => $create_extra_params{kilobytes_requested},
            );
        }
    } else {
        $allocation = $class->create(%params, %create_extra_params);
    }

    Genome::Sys->unlock_resource(resource_lock => $lock);

    unless ($allocation) {
        croak sprintf("Failed to shadow_get_or_create allocation from params %s",
                Data::Dumper::Dumper({%params, %create_extra_params}));
    }

    return $allocation;
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
        $params{id} = $class->id;
        $class = ref($class);
    }
    confess "Require allocation ID!" unless exists $params{id};

    my $allocation;
    if ($ENV{UR_DBI_NO_COMMIT}) {
        $allocation = $class->$method(%params);
    }
    else {
        # remove the parens if you DARE
        my @includes = map { ( '-I' => $_ ) } UR::Util::used_libs;

        my $param_string = Genome::Utility::Text::hash_to_string(\%params, 'q');
        my @statements = (
            qq(Genome::Utility::Instrumentation::timer('disk.allocation.require', sub { require $class })),
            sprintf('%s->%s(%s)', $class, $method, $param_string),
            q(UR::Context->commit),
        );
        my $perl_program_string = join('; ', @statements);

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
        $allocation = $class->_reload_allocation($params{id});
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

    my $allocation;
    Genome::Utility::Instrumentation::timer('disk.allocation.get_parent_allocation', sub {
        $allocation = $class->_get_parent_allocation_impl($path);
    });

    return $allocation if $allocation;
    return;
}

sub _get_parent_allocation_impl {
    my ($class, $path) = @_;
    Carp::confess("no path defined") unless defined $path;

    my ($allocation) = $class->get(allocation_path => $path);
    return $allocation if $allocation;

    my $dir = File::Basename::dirname($path);
    if ($dir ne '.' and $dir ne '/') {
        return $class->_get_parent_allocation_impl($dir);
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
        return !($class->get_child_allocations($path));
    }

    my ($err, $errstr, $row_arrayref);
    Genome::Utility::Instrumentation::timer('disk.allocation.child_allocation_query', sub {
        my $dbh = $data_source->get_default_handle();
        my $query_object = $dbh->prepare($query_string);
        $query_object->bind_param(1, $path . "/%");
        $query_object->execute();

        $row_arrayref = $query_object->fetchrow_arrayref();
        $err = $dbh->err;
        $errstr = $dbh->errstr;
        $query_object->finish();
    });

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

sub archive_after_time {
    my $self = shift;

    if (@_) {
        my ($new_time) = @_;
        my $old_time = $self->__archive_after_time();
        my $event = ($old_time lt $new_time) ? 'strengthened' : 'weakened';

        Genome::Timeline::Event::Allocation->$event(
            sprintf('original time: %s - new time: %s', $old_time, $new_time),
            $self,
        );
        return $self->__archive_after_time(@_);
    } else {
        return $self->__archive_after_time;
    }
}

sub purge {
    my $self = shift;
    my $reason = shift;
    die("You must supply a reason for the purge!") unless $reason;
    return $self->_execute_system_command('_purge', reason => $reason);
}

sub _purge {
    my $class = shift;
    my %params = @_;
    my $reason = delete $params{reason};
    my $allocation_id = delete $params{id};

    my $self = $class->get($allocation_id);
    die("No allocation found for id: $allocation_id") unless $self;
    die("You must supply a reason for the purge!") unless $reason;

    $self->_create_file_summaries();
    my $trash_folder = $self->_get_trash_folder();

    unless($ENV{UR_DBI_NO_COMMIT}) {
        my $destination_directory = Genome::Sys->create_directory(File::Spec->join($trash_folder, $self->id));
        dirmove($self->absolute_path, $destination_directory);
    }

    Genome::Timeline::Event::Allocation->purged(
        $reason,
        $self,
    );

    $self->status('purged');
    $self->kilobytes_requested(0);
    $self->kilobytes_used(0);

    return 1;
}

sub _create_file_summaries {
    my $self = shift;

    my $old_cwd = getcwd;
    chdir($self->absolute_path);
    my @files;
    #why is File::Find this stupid? who knows...
    find(sub { push(@files, $File::Find::name) unless (-d $_) }, '.');
    chdir($old_cwd);

    for my $file (@files) {
        Genome::Disk::Allocation::FileSummary->create_or_update(
            allocation => $self,
            file => $file
        );
    }
}

sub _symlink_new_path_from_old {
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

sub _default_archive_after_time {
    DateTime->now(time_zone => 'local')->add(years => 1)->strftime('%F %T');
}

sub _get_trash_folder {
    my $self = shift;

    my @dv = Genome::Disk::Volume->get(disk_group_names => 'apipe_trash');
    my %trash_map = map {
       $self->_extract_aggr($_->physical_path) => File::Spec->join($_->mount_path, '.trash');
    } @dv;

    my $aggr = $self->_extract_aggr($self->volume->physical_path);

    return $trash_map{$aggr};
}

sub _extract_aggr {
    my $self = shift;
    return (shift =~ m!/(aggr\d{2})/!)[0];
}

1;
