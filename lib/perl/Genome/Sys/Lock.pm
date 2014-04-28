package Genome::Sys::Lock;

use strict;
use warnings;

use Carp qw(carp croak);
use Genome::Sys::FileLock;
use Genome::Sys::Lock::NessyBackend;
use List::MoreUtils qw(all);

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
    my $class = shift;
    my %args = with_default_lock_resource_args(@_);

    my @locks;
    my $unwind = sub {
        for my $pair (@locks) {
            my ($backend, $resource_lock) = @$pair;
            my @unlock_args = $backend->translate_unlock_args(resource_lock => $resource_lock);
            $backend->unlock(@unlock_args);
        }
    };

    for my $backend (backends()) {
        my @lock_args = $backend->translate_lock_args(%args);
        my $lock = $backend->lock(@lock_args);
        if ($lock) {
            push @locks, [$backend, $lock];
        }

        if (is_mandatory($backend) && !$lock) {
            $unwind->();
            return;
        }
    }

    if (all { $_->has_lock($args{resource_lock}) } backends()) {
        Genome::Utility::Instrumentation::increment('genome.sys.lock.lock_resource.consistent');
    } else {
        Genome::Utility::Instrumentation::increment('genome.sys.lock.lock_resource.inconsistent');
    }

    $class->_cleanup_handler_check();

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
    my $class = shift;
    my %args = @_;

    my $rv = 1;
    for my $backend (backends()) {
        my @unlock_args = $backend->translate_unlock_args(%args);
        my $unlocked = $backend->unlock(@unlock_args);
        if (is_mandatory($backend)) {
            $rv = $rv && $unlocked;
        }
    }

    return $rv;
}

=item release_all()

C<release_all()> should release all locks managed by this process.  This should
not be called under normal circumstances.  Instead unlock each lock
individually.

=cut

sub release_all {
    for my $backend (backends()) {
        $backend->release_all();
    }
}

sub with_default_lock_resource_args {
    my %args = @_;

    $args{block_sleep} = 60 unless defined $args{block_sleep};
    $args{max_try} = 7200 unless defined $args{max_try};
    $args{wait_announce_interval} = 0 unless defined $args{wait_announce_interval};

    return %args;
}

my @backends = ('Genome::Sys::FileLock');
if ($ENV{GENOME_NESSY_SERVER}) {
    my $nessylock = Genome::Sys::Lock::NessyBackend->new(
        url => 'http://nessy.gsc.wustl.edu/',
        is_mandatory => 0,
    );
    push @backends, $nessylock;

    UR::Context->process->add_observer(
        aspect => 'sync_databases',
        callback => sub {
            my ($ctx, $aspect, $sync_db_result) = @_;
            if ($sync_db_result) {
                use vars '@CARP_NOT';
                local @CARP_NOT = (@CARP_NOT, 'UR::Context');
                foreach my $claim ($nessylock->claims) {
                    $claim->validate
                        || Carp::croak(sprintf('Claim %s failed to verify during commit', $claim->resource_name));
                }
            }
        }
    );
}

sub is_mandatory {
    my $backend = shift;
    return $backend->is_mandatory();
}

sub backends {
    return @backends;
}

sub add_backend {
    my ($class, $backend) = @_;
    push @backends, $backend;
}

sub remove_backend {
    my ($class, $backend) = @_;
    @backends = grep { $_ ne $backend } @backends;
}

########################################################################
# Private
########################################################################

my $_cleanup_handler_installed;
sub _cleanup_handler_check {
    my $class = shift;
    return if $_cleanup_handler_installed++;
    $SIG{'INT'} = \&_INT_cleanup;
    $SIG{'TERM'} = \&_INT_cleanup;
    $SIG{'HUP'} = \&_INT_cleanup;
    $SIG{'ABRT'} = \&_INT_cleanup;
    $SIG{'QUIT'} = \&_INT_cleanup;
    $SIG{'SEGV'} = \&_INT_cleanup;
}

sub _INT_cleanup {
    release_all();
    Carp::confess("INT/TERM cleanup activated in Genome::Sys::Lock");
}

END {
    release_all();
};

1;
