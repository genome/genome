package Genome::Sys::Lock;

use strict;
use warnings;

use Carp qw(carp croak);

use Genome::Sys;
use Genome::Sys::FileLock;
use Genome::Sys::NessyLock;

# backend => mandatory
my %backends = (
    'Genome::Sys::NessyLock' => 0,
    'Genome::Sys::FileLock' => 1
);

=item lock_resource()

Keyword Arguments:

=over

=item resource_lock

C<resource_lock> is the path to the file representing the lock for some
resource.

=item block_sleep

C<block_sleep> specifies the the number of seconds to sleep between attempts to
acquire the lock.

Default value: 60

=item max_try

C<max_try> specifies the number of retries to attempt to acquire the lock.

Default value: 7200

=item wait_announce_interval

C<wait_announce_interval> specifies the interval in seconds between status
updates indicating the contentious lock.  A zero value means to announce every attempt.

Default value: 0

=back

Returns:

=over

=item On Success

When C<lock()> acquired a lock on the resource it returns the lock path, i.e.
C<resource_lock>.

=item On Timeout

When C<lock()> is not able to acquire a lock on the resource within C<max_try>
(+1) attempts it returns C<undef>.

=item On Failure

When C<lock()> fails to create the lock, for reasons other than contention,
then C<lock()> will C<croak()>.

=back

=cut

sub lock_resource {
    my $self = shift;
    my %args = with_default_lock_resource_args(@_);

    my @locks;
    for my $backend (keys %backends) {
        my $mandatory = $backends{$backend};
        my $lock = $backend->lock(%args);
        if ($lock) {
            push @locks, $lock;
        }
        if ($mandatory && !$lock) {
            while (pop @locks) {
                $_->unlock(%args);
            }
            return;
        }
    }

    my $nessy_claim = Genome::Sys::NessyLock->has_lock($args{resource_lock});
    my $rv = Genome::Sys::FileLock->has_lock($args{resource_lock});

    if (Genome::Sys::NessyLock->is_enabled) {
        $self->_lock_resource_report_inconsistent_locks($args{resource_lock}, $rv, $nessy_claim);
    }

    $self->cleanup_handler_check();

    return $args{resource_lock};

}

=item unlock_resource()

Keyword Arguments:

=over

=item resource_lock

C<resource_lock> is the path to the file representing the lock for some
resource.

=back

=cut

sub unlock_resource {
    my $self = shift;
    my %args = @_;

    my $rv;
    for my $backend (keys %backends) {
        $rv = $backend->unlock(%args);
    }

    return $rv;
}

=item clear_state()

C<clear_state()> can be used after fork() to get a "clean" lock state.

=cut

sub clear_state {
    for my $backend (keys %backends) {
        $backend->clear_state();
    }
}

=item release_all()

C<release_all()> should release all locks managed by this process.  This should
not be called under normal circumstances.  Instead unlock each lock
individually.

=cut

sub release_all {
    for my $backend (keys %backends) {
        $backend->release_all();
    }
}

sub _lock_resource_report_inconsistent_locks {
    my($self, $resource_lock, $file_lock, $nessy_claim) = @_;

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
    my %args = @_;

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

sub with_default_lock_resource_args {
    my %args = @_;

    $args{block_sleep} = 60 unless defined $args{block_sleep};
    $args{max_try} = 7200 unless defined $args{max_try};
    $args{wait_announce_interval} = 0 unless defined $args{wait_announce_interval};

    return %args;
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

sub cleanup_handler_check {
    my $self = shift;
    $SIG{'INT'} = \&INT_cleanup;
    $SIG{'TERM'} = \&INT_cleanup;
    $SIG{'HUP'} = \&INT_cleanup;
    $SIG{'ABRT'} = \&INT_cleanup;
    $SIG{'QUIT'} = \&INT_cleanup;
    $SIG{'SEGV'} = \&INT_cleanup;
}

sub INT_cleanup {
    release_all();
    print STDERR "INT/TERM cleanup activated in Genome::Sys::Lock\n";
    Carp::confess;
}

END {
    release_all();
};

1;
