package Genome::Sys::Lock;

use strict;
use warnings;

use Carp qw(carp croak);
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

my %RESOURCE_LOCK_SCOPE;
sub lock_resource {
    my $class = shift;
    my %args = with_default_lock_resource_args(@_);
    validate_scope($args{scope});

    if (exists $RESOURCE_LOCK_SCOPE{$args{resource_lock}}
        && $RESOURCE_LOCK_SCOPE{$args{resource_lock}} ne $args{scope}) {
        Carp::confess(sprintf("Attempted to lock the resource (%s) in "
                . "scope (%s) while it has an existing lock in scope (%s).  "
                . "Locking in multiple scopes is not currently supported.",
                $args{resource_lock},
                $RESOURCE_LOCK_SCOPE{$args{resource_lock}},
                $args{scope}));
    }

    my @locks;
    my $unwind = sub {
        for my $pair (@locks) {
            my ($backend, $resource_lock) = @$pair;
            my @unlock_args = $backend->translate_unlock_args(resource_lock => $resource_lock);
            $backend->unlock(@unlock_args);
        }
    };

    for my $backend (backends($args{scope})) {
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

    if (all { $_->has_lock($args{resource_lock}) } backends($args{scope})) {
        Genome::Utility::Instrumentation::increment('genome.sys.lock.lock_resource.consistent');
    } else {
        Genome::Utility::Instrumentation::increment('genome.sys.lock.lock_resource.inconsistent');
    }

    $class->_cleanup_handler_check();

    $RESOURCE_LOCK_SCOPE{$args{resource_lock}} = $args{scope};

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

    my $scope = $RESOURCE_LOCK_SCOPE{$args{resource_lock}};

    my $rv = 1;
    for my $backend (backends($scope)) {
        if ($backend->has_lock($args{resource_lock})) {
            my @unlock_args = $backend->translate_unlock_args(%args);
            my $unlocked = $backend->unlock(@unlock_args);
            if (is_mandatory($backend)) {
                $rv = $rv && $unlocked;
            }
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
    for my $scope (scopes()) {
        for my $backend(backends($scope)) {
            $backend->release_all();
        }
    }
}

sub with_default_lock_resource_args {
    my %args = @_;

    $args{block_sleep} = 60 unless defined $args{block_sleep};
    $args{max_try} = 7200 unless defined $args{max_try};
    $args{wait_announce_interval} = 0 unless defined $args{wait_announce_interval};

    return %args;
}

sub is_mandatory {
    my $backend = shift;
    return $backend->is_mandatory();
}

my $backends = {};
sub backends {
    my $scope = shift;

    validate_scope($scope);

    if (exists $backends->{$scope}) {
        return @{$backends->{$scope}};
    } else {
        Carp::confess(sprintf("No backends registered for scope '%s'",
                $scope));
    }
}

sub all_backends {
    return %$backends;
}

sub scopes {
    return ('site', 'tgisan', 'unknown');
}

sub clear_backends {
    $backends = {};
}

sub set_backends {
    my $class = shift;
    my %new_backends = @_;

    $class->clear_backends();
    for my $scope (keys %new_backends) {
        for my $backend (@{$new_backends{$scope}}) {
            $class->add_backend($scope, $backend);
        }
    }
}

sub add_backend {
    my ($class, $scope, $backend) = @_;
    validate_scope($scope);

    push @{$backends->{$scope}}, $backend;
}

sub validate_scope {
    my $scope = shift;
    unless (grep {$scope eq $_} scopes()) {
        Carp::confess(sprintf("Invalid scope '%s'.  Valid scopes are [%s].",
                    $scope, join(', ', scopes())));
    }
}

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
