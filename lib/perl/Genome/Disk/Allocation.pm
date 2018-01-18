package Genome::Disk::Allocation;

use strict;
use warnings;

use Genome;
use Genome::Sys::LockProxy qw();
use Genome::Utility::Instrumentation;
use Genome::Utility::File::Mode qw(mode);

use Carp qw(croak confess);
use Digest::MD5 qw(md5_hex);
use File::Find;
use File::Find::Rule;
use Cwd;
use DateTime;
use Sub::Install;

require File::Spec;

our $TESTING_DISK_ALLOCATION = 0;

use constant SECONDS_IN_ONE_YEAR => 365*24*60*60;

class Genome::Disk::Allocation {
    roles => 'Genome::Role::Notable',
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
            len => 40,
            default_value => 'active',
            valid_values => [ "active", "purged", "archived", "invalid" ],
            doc => 'The current status of this allocation',
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
            doc => 'After this time, this allocation is subject to being archived',
        },
        timeline_events => {
            is => 'Genome::Timeline::Event::Allocation',
            reverse_as => 'object',
            where => [ -order_by => [ "created_at" ] ],
            is_many => 1,
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

sub set_read_only { shift->set_permissions_read_only(@_) }
sub set_permissions_read_only {
    my $self = shift;
    my @paths = File::Find::Rule->not(File::Find::Rule->symlink)->in($self->absolute_path);
    for my $path (@paths) {
        mode($path)->rm_all_writable();
    }
}

sub set_files_read_only {
    my $self = shift;
    my @paths = File::Find::Rule->not(File::Find::Rule->symlink)->file->in($self->absolute_path);
    for my $path (@paths) {
        mode($path)->rm_all_writable();
    }
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
        if ($CREATE_DUMMY_VOLUMES_FOR_TESTING) {
            my $tmp_volume = Genome::Disk::Volume->create_dummy_volume(
                mount_path => $params{mount_path},
                disk_group_name => $params{disk_group_name},
            );
            $params{mount_path} = $tmp_volume->mount_path;
        }
    }

    if ( Genome::Disk::Group->is_archive($params{disk_group_name}) ) {
        die $class->error_message('Cannot create disk allocation in an archive group: '.$params{disk_group_name});
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

sub copy {
    my ($class, %params) = @_;

    Genome::Utility::Instrumentation::inc('disk.allocation.copy');

    return $class->_execute_system_command('_copy', %params);
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

my $statuses = Genome::Disk::Allocation->__meta__->property_meta_for_name('status')->valid_values;
for my $status (@$statuses) {
    Sub::Install::install_sub({
        into => __PACKAGE__,
        as => "is_$status",
        code => sub { return shift->status eq $status; }
    });
}

sub archive_path {
    my $self = shift;
    my $mount_path = $self->volume->is_archive
                   ? $self->volume->mount_path
                   : $self->volume->archive_mount_path;
    return File::Spec->join(
        $mount_path,
        $self->group_subdirectory,
        $self->allocation_path,
    );
}

sub tar_path {
    my $self = shift;
    return File::Spec->join($self->archive_path, 'archive.tar');
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

    my %parameters = @_;
    $parameters{allocation_id} = delete $parameters{id};

    my $pars = Genome::Disk::Detail::Allocation::CreationParameters->create(
        %parameters);

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

    $self->_delete_timeline_events;
    my @deletion_observers = $self->_get_deletion_observers;
    $self->SUPER::delete;

    $class->_create_observer(@deletion_observers);

    return 1;
}

sub _delete_timeline_events {
    my $self = shift;

    for my $event ($self->timeline_events) {
        $event->delete();
    }
}

sub _get_deletion_observers {
    my $self = shift;

    my $class = $self->class;
    my $path = $self->absolute_path;

    if ($self->is_archived) {
        return sub {$class->_cleanup_archive_directory($path)};

    } else {
        return (
            $class->_mark_for_deletion_closure($path),
            $class->_remove_directory_closure($path),
        );
    }

}

sub _reallocate {
    my $class = shift;

    my %parameters = @_;
    $parameters{allocation_id} = delete $parameters{id};

    my $reallocator = Genome::Disk::Detail::Allocation::Reallocator->create(
        %parameters);

    return $reallocator->reallocate;
}

sub _move {
    my $class = shift;

    my %parameters = @_;
    $parameters{allocation_id} = delete $parameters{id};

    my $mover = Genome::Disk::Detail::Allocation::Mover->create(
        %parameters);

    return $mover->move;
}

sub _copy {
    my $class = shift;

    my %parameters = @_;
    $parameters{allocation_id} = delete $parameters{id};

    my $copier = Genome::Disk::Detail::Allocation::Copier->create(
        %parameters);

    return $copier->copy;
}

sub _archive {
    my $class = shift;

    my %parameters = @_;
    $parameters{allocation_id} = delete $parameters{id};

    my $archiver = Genome::Disk::Detail::Allocation::Archiver->create(
        %parameters);

    return $archiver->archive;
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
    my $class = shift;

    my %parameters = @_;
    $parameters{allocation_id} = delete $parameters{id};

    my $unarchiver = Genome::Disk::Detail::Allocation::Unarchiver->create(
        %parameters);

    return $unarchiver->unarchive;
}

# Locks the allocation, if lock is not manually released (it had better be!)
# it'll be automatically cleaned up on program exit
sub get_lock {
    my ($class, $id, $tries) = @_;
    $tries ||= 60;

    my $allocation_lock = Genome::Sys::LockProxy->new(
        resource => 'allocation/allocation_' . join('_', split(' ', $id)),
        scope => 'site',
    )->lock(
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
        $lock->unlock();
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

    my $resource = File::Spec->join(
        'allocation',
        'allocation_path_' . $md5_hex,
    );

    my $lock = Genome::Sys::LockProxy->new(
        resource => $resource,
        scope => 'site',
    )->lock(
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

    $lock->unlock();

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

sub __rollback_property__ {
    my ($self, $property_name) = @_;

    # We do not want to create new Genome::Timeline::Events during rollback so
    # have to bypass the overridden methods and directly use the accessors.
    if ($property_name eq 'archivable') {
        my $saved_value = UR::Context->current->value_for_object_property_in_underlying_context($self, $property_name);
        return $self->__archivable($saved_value);
    }
    if ($property_name eq 'archive_after_time') {
        my $saved_value = UR::Context->current->value_for_object_property_in_underlying_context($self, $property_name);
        return $self->__archive_after_time($saved_value);
    }

    return $self->SUPER::__rollback_property__($property_name);
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
        my @includes = map { -I => $_ } UR::Util::used_libs;

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
        if ($allocation) {
            my $owner = $allocation->owner;
            if($owner and UR::Context->current->object_exists_in_underlying_context($owner)) {
                UR::Context->current->reload($owner);
            }
        }
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

    Genome::Sys::CommitAction->create( on_commit => $callback );
}

# Returns a closure that creates a directory at the given path
sub _create_directory_closure {
    my ($class, $path) = @_;
    return sub {
        # This method currently returns the path if it already exists instead of failing
        my $dir = eval { Genome::Sys->create_directory($path) };
        unless (defined $dir and -d $dir) {
            $class->error_message("Could not create allocation directory at $path!\n$@");
        }
    };
}

# Returns a closure that removes the given directory
sub _remove_directory_closure {
    my ($class, $path) = @_;
    return sub {
        if (-d $path and not $ENV{UR_DBI_NO_COMMIT}) {
            $class->debug_message("Removing allocation directory $path");
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
            $class->debug_message("Marking directory at $path as deallocated");
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
            mode($File::Find::name)->rm_all_writable();
        };

        $class->debug_message("Marking directory at $path read-only");
        File::Find::find(\&mark_read_only, $path);
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

    return 1 if $ENV{UR_DBI_NO_COMMIT} and not $TESTING_DISK_ALLOCATION;

    $path =~ s/\/+$//;

    my $meta        = $class->__meta__;
    my $table_name  = $meta->table_name;
    my $data_source = $meta->data_source;
    my $owner       = $data_source->owner;
    my $query_string;

    if ($data_source->isa('UR::DataSource::Oracle')) {
        my $fq_table_name = join('.', $owner, $table_name);
        $query_string = sprintf(
            q(select 1 from %s where allocation_path like ? AND rownum <= 1),
            $fq_table_name);
    } elsif ($data_source->isa('UR::DataSource::Pg') || $data_source->isa('UR::DataSource::SQLite')) {
        $query_string = sprintf(
            q(select 1 from %s where allocation_path like ? LIMIT 1),
            $table_name);
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
    my $cmd = "if [ -d $directory ] ; then rm -rf $directory ; else exit 1; fi";
    unless ($ENV{UR_DBI_NO_COMMIT}) {
        my $guard = Genome::Config::set_env('lsb_sub_additional', ''); #no docker for archives
        my ($job_id, $status) = Genome::Sys->bsub_and_wait(
            queue => Genome::Config::get('archive_lsf_queue'),
            cmd => "$cmd",
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

sub _resolve_param_value_from_text_by_name_or_id {
    my ($class, $param_arg) = @_;

    my @results = Command::V2->_resolve_param_value_from_text_by_name_or_id($class, $param_arg);

    if(!@results) {
        my $allocation = $class->get_allocation_for_path($param_arg);
        push @results, $allocation if $allocation;
    }

    return @results;
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

    my %parameters = @_;
    $parameters{allocation_id} = delete $parameters{id};

    my $purger = Genome::Disk::Detail::Allocation::Purger->create(
        %parameters);

    return $purger->purge;
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
    DateTime->now(time_zone => 'local')->add(seconds => SECONDS_IN_ONE_YEAR)->strftime('%F 00:00:00');
}

1;
