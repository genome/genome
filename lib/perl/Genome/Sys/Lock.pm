package Genome::Sys::Lock;
use strict;
use warnings;
use Time::HiRes;
use File::Basename;
use MIME::Lite;
use Carp qw(carp croak);
use File::Temp;
use Sys::Hostname qw/hostname/;
use Genome;
use base 'UR::ModuleBase';   # *_message methods, but no constructor
use Genome::Utility::Instrumentation;

my %SYMLINKS_TO_REMOVE;
my %NESSY_LOCKS_TO_REMOVE;
my $LOCKING_CLIENT;

sub lock_resource {
    my ($self,%args) = @_;

    $args{block_sleep} = 60 unless defined $args{block_sleep};
    $args{max_try} = 7200 unless defined $args{max_try};
    $args{wait_announce_interval} = 0 unless defined $args{wait_announce_interval};

    @args{'resource_lock', 'parent_dir'} = $self->_resolve_resource_lock_and_parent_dir_for_lock_resource(%args);

    $self->_start_locking_client;

    my $nessy_claim = $self->_new_style_lock(%args);
    my $rv = $self->_file_based_lock_resource(%args);
    $self->_lock_resource_report_inconsistent_locks($args{resource_lock}, $rv, $nessy_claim);

    return $rv;

}

sub _lock_resource_report_inconsistent_locks {
    my($self, $resource_lock, $file_lock, $nessy_claim) = @_;

    return unless $LOCKING_CLIENT;

    my $t = "%s-lock acquired but %s-based did not: $resource_lock";

    if ($nessy_claim and !$file_lock) {
        Genome::Utility::Instrumentation::increment('sys.lock.locked_nessy_only');
        Genome::Logger->debugf($t, 'Nessy', 'File');
        return;
    }

    if ($file_lock and !$nessy_claim) {
        Genome::Utility::Instrumentation::increment('sys.lock.locked_file_only');
        Genome::Logger->debugf($t, 'File', 'Nessy');
        return;
    }

    if ($file_lock and $nessy_claim) {
        Genome::Utility::Instrumentation::increment('sys.lock.locked_both');
        return;
    }
}


sub _resolve_resource_lock_and_parent_dir_for_lock_resource {
    my($self, %args) = @_;

    my $resource_lock = delete $args{resource_lock};
    my ($lock_directory,$resource_id,$parent_dir);
    if ($resource_lock) {
        $parent_dir = File::Basename::dirname($resource_lock);
        Genome::Sys->create_directory($parent_dir);
        unless (-d $parent_dir) {
            Carp::croak("failed to make parent directory $parent_dir for lock $resource_lock!: $!");
        }
    }
    else {
        $lock_directory =  delete $args{lock_directory} || Carp::croak('Must supply lock_directory to lock resource');
        Genome::Sys->create_directory($lock_directory);
        $resource_id = $args{'resource_id'} || Carp::croak('Must supply resource_id to lock resource');
        $resource_lock = $lock_directory . '/' . $resource_id . ".lock";
        $parent_dir = $lock_directory
    }

    return ($resource_lock, $parent_dir);
}

sub _file_based_lock_resource {
    my($self, %args) = @_;

    my $total_lock_start_time = Time::HiRes::time();

    my($resource_lock, $parent_dir)
        = delete @args{'resource_lock','parent_dir'};

    my $basename = File::Basename::basename($resource_lock);

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

    my $owner_details = $self->_resolve_lock_owner_details;
    my $lock_dir_template = sprintf("lock-%s--%s_XXXX",$basename,$owner_details);
    my $tempdir =  File::Temp::tempdir($lock_dir_template, DIR=>$parent_dir, CLEANUP=>1);

    unless (-d $tempdir) {
        Carp::croak("Failed to create temp lock directory ($tempdir)");
    }

    # make this readable for everyone
    chmod(0777, $tempdir) or Carp::croak("Can't chmod 0777 path ($tempdir): $!");

    # drop an info file into here for compatibility's sake with old stuff.
    # put a "NOKILL" here on LSF_JOB_ID so an old process doesn't try to snap off the job ID and kill me.
    my $lock_info = IO::File->new("$tempdir/info", ">");
    unless ($lock_info) {
        Carp::croak("Can't create info file $tempdir/info: $!");
    }

    my ($my_host, $my_pid, $my_lsf_id, $my_user) = (hostname, $$, ($ENV{'LSB_JOBID'} || 'NONE'), Genome::Sys->username);
    my $job_id = (defined $ENV{'LSB_JOBID'} ? $ENV{'LSB_JOBID'} : "NONE");
    $lock_info->printf("HOST %s\nPID $$\nLSF_JOB_ID_NOKILL %s\nUSER %s\n",
                       $my_host,
                       $ENV{'LSB_JOBID'},
                       $ENV{'USER'},
                     );
    $lock_info->close();

    my $initial_time = time;
    my $last_wait_announce_time = $initial_time;

    my $lock_attempts = 1;
    my $ret;
    while(!($ret = symlink($tempdir,$resource_lock))) {
        # TONY: The only allowable failure is EEXIST, right?
        # If any other error comes through, we end up in bigger trouble.
        use Errno qw(EEXIST ENOENT :POSIX);
        if ($! != EEXIST) {
            $self->error_message("Can't create symlink from $tempdir to lock resource $resource_lock because: $!");
            Carp::croak($self->error_message());
        }
        my $symlink_error = $!;
        chomp $symlink_error;
        return undef unless $max_try--;

        my $target = readlink($resource_lock);

        #If symlink no longer exists
        if ($! == ENOENT || !$target) {
            sleep $block_sleep;
            next;
        }

        #symlink could have disappeared between the top of the while loop and now, so we check to see if it's changed
        my $target_exists = (-e $target);
        unless($target eq readlink($resource_lock)){
            next;
        }

        if (!$target_exists) {
            # TONY: This means the lock symlink points to something that's been deleted
            # That's _really_ bad news and should probably get an email like below.
            $self->error_message("Lock ($resource_lock) exists but target ($target) does not exist.");
            Carp::croak($self->error_message);
        }

        if ($self->is_my_lock_target($target)) {
            $self->error_message("Tried to lock resource more than once: $resource_lock");
            Carp::croak($self->error_message);
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
            $self->status_message("waiting (total_elapsed_time = $total_elapsed_time seconds) on lock for resource '$resource_lock': $symlink_error. lock_info is:\n$info_content");
        }

        if ($lsf_id ne "NONE") {
            my ($job_info,$events) = Genome::Model::Event->lsf_state($lsf_id);
                 unless ($job_info) {
                     Genome::Utility::Instrumentation::increment('sys.lock.lock.found_orphan');
                     $self->warning_message("Invalid lock for resource $resource_lock\n"
                                            ." lock info was:\n". $info_content ."\n"
                                            ."Removing old resource lock $resource_lock\n");
                     unless ($Genome::Sys::IS_TESTING) {
                        my $message_content = <<END_CONTENT;
Hey Apipe,

This is a lock attempt on %s running as PID $$ LSF job %s and user %s.

I'm about to remove a lock file held by a LSF job that I think is dead.

The resource is:

%s

Here's info about the job that I think is gone.

%s

I'll remove the lock in an hour.  If you want to save the lock, kill me
before I unlock the process!

Your pal,
Genome::Utility::Filesystem

END_CONTENT

                        my $msg = MIME::Lite->new(From    => sprintf('"Genome::Utility::Filesystem" <%s@%s>', $ENV{'USER'}, $ENV{GENOME_EMAIL_DOMAIN}),
                                              To      => $ENV{GENOME_EMAIL_PIPELINE_NOISY},
                                              Subject => 'Attempt to release a lock held by a dead process',
                                              Data    => sprintf($message_content, $my_host, $ENV{'LSB_JOBID'}, $ENV{'USER'}, $resource_lock, $info_content),
                                            );
                        $msg->send();
                        $self->status_message('Sleeping for one hour...');
                        sleep 60 * 60;
                 }
                     $self->unlock_resource(resource_lock => $resource_lock, force => 1);
                     #maybe warn here before stealing the lock...
               }
           }
        sleep $block_sleep;
        $lock_attempts += 1;
       }
    $SYMLINKS_TO_REMOVE{$resource_lock} = 1;

    # do we need to activate a cleanup handler?
    $self->cleanup_handler_check();

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

sub _start_locking_client {
    my $class = shift;

    if ($ENV{GENOME_NESSY_SERVER} and ! $LOCKING_CLIENT) {
        require Nessy::Client;
        $LOCKING_CLIENT = Nessy::Client->new(url => $ENV{GENOME_NESSY_SERVER});
    }
}

sub _new_style_lock {
    my($self, %args) = @_;

    return 1 unless $LOCKING_CLIENT;

    my %user_data;
    @user_data{'host','pid','lsf_id','user'}
        = (hostname, $$, ($ENV{'LSB_JOBID'} || 'NONE'), Genome::Sys->username);

    my $resource_lock = $args{resource_lock};

    if ($self->_is_holding_nessy_lock($resource_lock)) {
        $self->error_message("Tried to lock resource more than once: $resource_lock");
        Carp::croak($self->error_message);
    }

    my $timeout = $self->_new_style_lock_timeout_from_args(%args);
    my $wait_announce_interval = delete $args{wait_announce_interval};
    unless (defined $wait_announce_interval) {
        croak('wait_announce_interval not defined');
    }


    my $claim = $LOCKING_CLIENT->claim($resource_lock, timeout => $timeout, user_data => \%user_data);
    $NESSY_LOCKS_TO_REMOVE{$resource_lock} = $claim if $claim;
    return $claim;
}

sub min_timeout {
    return 5;
}
sub _new_style_lock_timeout_from_args {
    my($self, %args) = @_;

    my $block_sleep = delete $args{block_sleep} || 0;

    my $max_try = delete $args{max_try} || 0;

    my $min_timeout = min_timeout();
    my $timeout = $max_try * $block_sleep;
    unless ($timeout > $min_timeout) {
        $timeout = $min_timeout;
        carp("increasing timeout to minimum ($min_timeout)");
    }

    return $timeout;
}

sub _is_holding_nessy_lock {
    my($self, $resource_lock) = @_;
    return $NESSY_LOCKS_TO_REMOVE{$resource_lock};
}

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

sub is_my_lock_target {
    my $self = shift;
    my $target = shift;

     my $basename = File::Basename::basename($target);
    $basename =~ s/_.{4}$//; #remove random chars tmpdir adds to end
    my $expected_details = $self->_resolve_lock_owner_details;
    my $length = length($expected_details);

    return (substr($basename, (0-$length)) eq $expected_details);
}


sub unlock_resource {
    my ($self,%args) = @_;

    my $resource_lock = $self->_resolve_resource_lock_for_unlock_resource(%args);

    my $rv = $self->_file_based_unlock_resource(
            resource_lock => $resource_lock,
            %args,
        );

    unless ($self->_new_style_release($resource_lock)) {
        $self->error_message("file-based released the lock, but Nessy did not.  resource_lock: $resource_lock");
    }
    return $rv;
}

sub _resolve_resource_lock_for_unlock_resource {
    my($self, %args) = @_;

    my $resource_lock = $args{resource_lock};
    unless ($resource_lock) {
        my ($lock_directory,$resource_id);
        $lock_directory =  delete $args{lock_directory} || Carp::croak('Must supply lock_directory to lock resource');
        $resource_id = $args{'resource_id'} || Carp::croak('Must supply resource_id to lock resource');
        $resource_lock = $lock_directory . '/' . $resource_id . ".lock";
    }
    return $resource_lock;
}

sub _file_based_unlock_resource {
    my($self, %args) = @_;

    my $resource_lock = delete $args{resource_lock};
    my $force = delete $args{force};

    my $target = readlink($resource_lock);
    if (!$target) {
        if ($! == ENOENT) {
            Genome::Utility::Instrumentation::increment('sys.lock.unlock.stolen_from_me');
            $self->error_message("Tried to unlock something that's not locked -- $resource_lock.");
            Carp::croak($self->error_message);
        } else {
            $self->error_message("Couldn't readlink $resource_lock: $!");
        }
    }
    unless (-d $target) {
        $self->error_message("Lock symlink '$resource_lock' points to something that's not a directory - $target. ");
        Carp::croak($self->error_message);
    }

    unless ($force) {
        unless ($self->is_my_lock_target($target)) {
             my $basename = File::Basename::basename($target);
             my $expected_details = $self->_resolve_lock_owner_details;
             Genome::Utility::Instrumentation::increment('sys.lock.unlock.stolen_from_me');
             $self->error_message("This lock does not look like it belongs to me.  $basename does not match $expected_details.");
             delete $SYMLINKS_TO_REMOVE{$resource_lock}; # otherwise the lock would be forcefully cleaned up when process exits
             Carp::croak($self->error_message);
        }
    }

    my $unlink_rv = unlink($resource_lock);
    if (!$unlink_rv) {
        Genome::Utility::Instrumentation::increment('sys.lock.unlock.unlink_failed');
        $self->error_message("Failed to remove lock symlink '$resource_lock':  $!");
        Carp::croak($self->error_message);
    }

    my $rmdir_rv = File::Path::rmtree($target);
    if (!$rmdir_rv) {
        Genome::Utility::Instrumentation::increment('sys.lock.unlock.rmtree_failed');
        $self->error_message("Failed to remove lock symlink target '$target', but we successfully unlocked.");
        Carp::croak($self->error_message);
    }

    delete $SYMLINKS_TO_REMOVE{$resource_lock};
    $self->cleanup_handler_check();

    Genome::Utility::Instrumentation::increment('sys.lock.unlock.success');
    return 1;
}

sub _new_style_release {
    my($self, $resource_lock) = @_;

    if ($LOCKING_CLIENT) {
        my $claim = delete $NESSY_LOCKS_TO_REMOVE{$resource_lock};
        if ($claim) {
            $claim->release;
        } else {
            $self->error_message("Nessy tried to release, but no claim in slot for resource_lock: $resource_lock");
        }
    } else {
        return 1;
    }
}

# FIXME - I think this is a private function to Filesystem.pm
sub cleanup_handler_check {
    my $self = shift;

    my $symlink_count = scalar keys %SYMLINKS_TO_REMOVE;

    if ($symlink_count > 0) {
        $SIG{'INT'} = \&INT_cleanup;
        $SIG{'TERM'} = \&INT_cleanup;
        $SIG{'HUP'} = \&INT_cleanup;
        $SIG{'ABRT'} = \&INT_cleanup;
        $SIG{'QUIT'} = \&INT_cleanup;
        $SIG{'SEGV'} = \&INT_cleanup;
    } else {
        delete $SIG{'INT'};
        delete $SIG{'TERM'};
        delete $SIG{'HUP'};
        delete $SIG{'ABRT'};
        delete $SIG{'QUIT'};
        delete $SIG{'SEGV'};
    }
}

END {
    exit_cleanup();
};

sub INT_cleanup {
    exit_cleanup();
    print STDERR "INT/TERM cleanup activated in Genome::Utility::Filesystem\n";
    Carp::confess;
}

sub exit_cleanup {
    for my $sym_to_remove (keys %SYMLINKS_TO_REMOVE) {
        if (-l $sym_to_remove) {
            warn("Removing remaining resource lock: '$sym_to_remove'") unless $ENV{'HARNESS_ACTIVE'};
            unlink($sym_to_remove) or warn "Can't unlink $sym_to_remove: $!";
        }
    }
    if ($LOCKING_CLIENT) {
        foreach my $resource_lock ( keys %NESSY_LOCKS_TO_REMOVE ) {
            __PACKAGE__->_new_style_release($resource_lock);
        }
        %NESSY_LOCKS_TO_REMOVE = ();
        undef $LOCKING_CLIENT;
    }
}

UR::Context->process->add_observer(
    aspect => 'sync_databases',
    callback => sub {
        my($ctx, $aspect, $sync_db_result) = @_;
        if ($sync_db_result) {
            use vars '@CARP_NOT';
            local @CARP_NOT = (@CARP_NOT, 'UR::Context');
            foreach my $claim (values %NESSY_LOCKS_TO_REMOVE ) {
                $claim->validate
                    || Carp::croak(sprintf('Claim %s failed to verify during commit', $claim->resource_name));
            }
        }
    }
);

1;
