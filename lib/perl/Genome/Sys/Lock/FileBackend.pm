package Genome::Sys::Lock::FileBackend;

use strict;
use warnings;

use Carp qw(carp croak);
use File::Basename;
use File::Temp;
use MIME::Lite;
use Sys::Hostname qw(hostname);
use Time::HiRes;
use Path::Class qw();

use Genome::Logger;
use Genome::Utility::Instrumentation;

use Mouse;
with qw(Genome::Sys::Lock::Backend);

has 'parent_dir' => (is => 'ro', isa => 'Str', required => 1);
has 'owned_resources' => (is => 'rw', isa => 'HashRef', lazy_build => 1);


sub lock {
    my($self, %args) = @_;

    my $total_lock_start_time = Time::HiRes::time();

    my $resource_lock = delete $args{'resource_lock'};
    unless ($resource_lock) {
        croak('resource_lock not set');
    }

    my $block_sleep = delete $args{block_sleep};
    unless (defined $block_sleep) {
        croak('block_sleep not defined');
    }

    my $max_try = delete $args{max_try};
    unless (defined $max_try) {
        croak('max_try not defined');
    }

    my $wait_announce_interval = delete $args{wait_announce_interval};
    unless (defined $wait_announce_interval) {
        croak('wait_announce_interval not defined');
    }

    my $symlink_path = $self->path_for_resource($resource_lock);

    my $tempdir = $self->tempdir_for_resource($resource_lock);

    $self->write_lock_info($tempdir);

    my $initial_time = time;
    my $last_wait_announce_time = $initial_time;

    my $lock_attempts = 1;
    my $ret;
    while(!($ret = symlink($tempdir,$symlink_path))) {
        # TONY: The only allowable failure is EEXIST, right?
        # If any other error comes through, we end up in bigger trouble.
        use Errno qw(EEXIST ENOENT :POSIX);
        if ($! != EEXIST) {
            Genome::Logger->fatal("Can't create symlink from $tempdir to lock resource $symlink_path because: $!\n");
        }
        my $symlink_error = $!;
        chomp $symlink_error;
        return undef unless $max_try--;

        my $target = readlink($symlink_path);

        #If symlink no longer exists
        if ($! == ENOENT || !$target) {
            sleep $block_sleep;
            next;
        }

        #symlink could have disappeared between the top of the while loop and now, so we check to see if it's changed
        my $target_exists = (-e $target);
        unless($target eq readlink($symlink_path)){
            next;
        }

        if (!$target_exists) {
            # TONY: This means the lock symlink points to something that's been deleted
            # That's _really_ bad news and should probably get an email like below.
            Genome::Logger->fatal("Lock ($resource_lock) exists but target ($target) does not exist.\n");
        }

        if ($self->is_my_lock_target($target)) {
            Genome::Logger->fatal("Tried to lock resource more than once: $resource_lock\n");
        }

        my $target_basename = File::Basename::basename($target);
        $target_basename =~ s/lock-.*?--//; #might not always work if lock name starts with "-" or contains "--"
        my ($host, $user, $pid, $lsf_id) = split /_/, $target_basename;
        my $info_content=sprintf("HOST %s\nPID %s\nLSF_JOB_ID %s\nUSER %s",$host,$pid,$lsf_id,$user);

        my $time = time;
        my $elapsed_time = $time - $last_wait_announce_time;
        if ($elapsed_time >= $wait_announce_interval) {
            $last_wait_announce_time = $time;
            my $total_elapsed_time = $time - $initial_time;
            Genome::Logger->notice("waiting (total_elapsed_time = $total_elapsed_time seconds) on lock for resource '$resource_lock': $symlink_error. lock_info is:\n$info_content\n");
        }

        if ($lsf_id ne "NONE") {
            my ($job_info,$events) = Genome::Model::Event->lsf_state($lsf_id);
            unless ($job_info) {
                Genome::Utility::Instrumentation::increment('sys.lock.lock.found_orphan');
                Genome::Logger->warning("Invalid lock for resource $resource_lock\n"
                    ." lock info was:\n". $info_content ."\n"
                    ."Removing old resource lock $resource_lock\n");
                $self->unlock(resource_lock => $resource_lock, force => 1);
            }
        }
        sleep $block_sleep;
        $lock_attempts += 1;
    }
    $self->owned_resources->{$resource_lock} = $$;

    my $total_lock_stop_time = Time::HiRes::time();
    my $lock_time_miliseconds = 1000 * ($total_lock_stop_time - $total_lock_start_time);

    my $caller_name = _resolve_caller_name(caller());

    Genome::Utility::Instrumentation::timing("lock_resource.$caller_name", $lock_time_miliseconds);
    Genome::Utility::Instrumentation::timing("lock_resource_attempts.$caller_name",
        $lock_attempts);

    Genome::Utility::Instrumentation::timing('lock_resource.total', $lock_time_miliseconds);
    Genome::Utility::Instrumentation::timing('lock_resource_attempts.total',
        $lock_attempts);

    # I don't think the conditional is actually needed but being safe
    if ($resource_lock) {
        Genome::Utility::Instrumentation::increment('sys.lock.lock.success');
    }

    return $resource_lock;
}

sub unlock {
    my($self, %args) = @_;

    my $resource_lock = delete $args{resource_lock};
    unless ($resource_lock) {
        carp('resource_lock is not set');
    }
    my $force = delete $args{force};

    my $symlink_path = $self->path_for_resource($resource_lock);

    my $target = readlink($symlink_path);
    if (!$target) {
        if ($! == ENOENT) {
            Genome::Utility::Instrumentation::increment('sys.lock.unlock.stolen_from_me');
            Genome::Logger->fatal("Tried to unlock something that's not locked -- $resource_lock.\n");
        } else {
            Genome::Logger->error("Couldn't readlink $symlink_path: $!\n");
        }
    }
    unless (-d $target) {
        Genome::Logger->fatal("Lock symlink '$symlink_path' points to something that's not a directory - $target. \n");
    }

    unless ($force) {
        unless ($self->is_my_lock_target($target)) {
             my $basename = File::Basename::basename($target);
             my $expected_details = $self->_resolve_lock_owner_details;
             Genome::Utility::Instrumentation::increment('sys.lock.unlock.stolen_from_me');
             delete $self->owned_resources->{$resource_lock}; # otherwise the lock would be forcefully cleaned up when process exits
             Genome::Logger->fatal("This lock does not look like it belongs to me.  $basename does not match $expected_details.\n");
        }
    }

    my $unlink_rv = unlink($symlink_path);
    if (!$unlink_rv) {
        Genome::Utility::Instrumentation::increment('sys.lock.unlock.unlink_failed');
        Genome::Logger->fatal("Failed to remove lock symlink '$symlink_path':  $!\n");
    }

    my $rmdir_rv = File::Path::rmtree($target);
    if (!$rmdir_rv) {
        Genome::Utility::Instrumentation::increment('sys.lock.unlock.rmtree_failed');
        Genome::Logger->fatal("Failed to remove lock symlink target '$target', but we successfully unlocked.\n");
    }

    delete $self->owned_resources->{$resource_lock};

    Genome::Utility::Instrumentation::increment('sys.lock.unlock.success');
    return 1;
}

sub clear_state {
    my $self = shift;
    $self->owned_resources = {};
}

sub translate_lock_args { shift; return @_ }
sub translate_unlock_args { shift; return @_ }

sub _resolve_caller_name {
    my ($package, $filename, $line) = @_;
    $package =~ s/::/./g;
    return $package;
}

#build the string for locks generated by this process
sub _resolve_lock_owner_details {
    my $self = shift;
    my ($my_host, $my_pid, $my_lsf_id, $my_user) = (hostname, $$, ($ENV{'LSB_JOBID'} || 'NONE'), Genome::Sys->username);
    my $job_id = (defined $ENV{'LSB_JOBID'} ? $ENV{'LSB_JOBID'} : "NONE");
    my $lock_details = join('_',$my_host,$ENV{'USER'},$$,$job_id);

    return $lock_details;
}

sub has_lock {
    my ($self, $resource_lock) = @_;
    return exists $self->owned_resources->{$resource_lock};
}

sub is_my_lock_target {
    my $self = shift;
    my $target = shift;

     my $basename = File::Basename::basename($target);
    $basename =~ s/_.{4}$//; #remove random chars tmpdir adds to end
    my $expected_details = $self->_resolve_lock_owner_details;
    my $length = length($expected_details);

    return (substr($basename, (0-$length)) eq $expected_details);
}

sub release_all {
    my $self = shift;

    my %errors;
    for my $owned_resource (keys %{$self->owned_resources}) {
        if ($self->owned_resources->{$owned_resource} == $$) {
            if (-l $self->path_for_resource($owned_resource)) {
                warn("Removing remaining resource lock: '$owned_resource'") unless $ENV{'HARNESS_ACTIVE'};
                eval {$self->unlock(resource_lock => $owned_resource)};
                $errors{$owned_resource} = $@ if $@;
            }
        } else {
            Genome::Logger->warning(sprintf(
                "Not cleaning up lock named (%s) because ".
                "the pid that created it (%s) is not my pid (%s).\n",
                $owned_resource,
                $self->owned_resources->{$owned_resource},
                $$,
            ));
        }
    }
    if (keys %errors) {
        Genome::Logger->error(sprintf(
            "Cleaning up the following locks failed: %s\n",
            join(", ", keys %errors),
        ));
        Carp::croak(join("\n", values %errors));
    }
}

sub path_for_resource {
    my ($self, $resource_lock) = @_;

    return Path::Class::file($self->parent_dir, $resource_lock)->stringify;
}

sub lock_dir_for_resource {
    my ($self, $resource_lock) = @_;

    return Path::Class::file($self->parent_dir, $resource_lock)->parent->stringify;
}

sub make_dir_path {
    my ($self, $dir_path) = @_;
    my $obj = Path::Class::dir($dir_path);
    $obj->mkpath(0, 0777);
}

sub tempdir_for_resource {
    my ($self, $resource_lock) = @_;

    my $lock_dir = $self->lock_dir_for_resource($resource_lock);
    $self->make_dir_path($lock_dir);

    my $basename = File::Basename::basename($resource_lock);

    my $owner_details = $self->_resolve_lock_owner_details;
    my $lock_dir_template = sprintf("lock-%s--%s_XXXX", $basename,
        $owner_details);
    my $tempdir =  File::Temp::tempdir($lock_dir_template,
        DIR => $lock_dir, CLEANUP => 1);

    unless (-d $tempdir) {
        Carp::croak("Failed to create temp lock directory ($tempdir)");
    }

    # make this readable for everyone
    chmod(0770, $tempdir) or Carp::croak("Can't chmod 0770 path ($tempdir): $!");

    return $tempdir;
}

sub write_lock_info {
    my ($self, $tempdir) = @_;

    # drop an info file into here for compatibility's sake with old stuff.
    # put a "NOKILL" here on LSF_JOB_ID so an old process doesn't try to snap off the job ID and kill me.
    my $lock_info = IO::File->new("$tempdir/info", ">");
    unless ($lock_info) {
        Carp::croak("Can't create info file $tempdir/info: $!");
    }

    my ($my_host, $my_pid, $my_lsf_id, $my_user) = (hostname, $$, ($ENV{'LSB_JOBID'} || 'NONE'), Genome::Sys->username);
    $lock_info->printf("HOST %s\nPID $$\nLSF_JOB_ID_NOKILL %s\nUSER %s\n",
                       $my_host,
                       $ENV{'LSB_JOBID'},
                       $ENV{'USER'},
                     );
    $lock_info->close();
}

sub _build_owned_resources {
    my $self = shift;
    return {};
}

__PACKAGE__->meta->make_immutable();
