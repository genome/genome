package Genome::Sys::LockProxy;

use strict;
use warnings;

use Carp qw(croak);
use Genome::Sys::Lock qw();
use Params::Validate qw(validate_with);
use Scalar::Util qw(blessed);

require Scope::Guard;

=item new()

Keyword Arguments:

=over

=item resource

C<resource> is a unique identifier for some resource.

=item scope

C<scope> is the scope to which the lock is bound.  See C<Genome::Sys::Lock::scopes()> for valid scopes.

=cut

sub new {
    my $class = shift;
    if (ref $class) {
        croak 'new() should be called as a class method';
    }
    my %params = validate_with(
        params => \@_,
        spec => Genome::Sys::Lock::PROXY_CONSTRUCTOR_PARAMS_SPEC(),
    );
    return bless \%params, $class;
}

=item resource

Accesor for the C<resource> field.

=cut

sub resource {
    my $self = shift;
    return $self->{resource};
}

=item scope

Accesor for the C<scope> field.

=cut

sub scope {
    my $self = shift;
    return $self->{scope};
}

=item lock()

C<lock()> is an instance method with optional parameters C<block_sleep>,
C<max_try>, or C<wait_announce_interval>.  It behaves similar to
C<Genome::Sys::Lock::lock_resource()> but on sucess C<lock()> returns an
instance of C<Genome::Sys::LockProxy> that can later be unlocked with
C<unlock()>.

=cut

sub lock {
    my $self = shift;
    my %params = validate_with(
        params => \@_,
        spec => Genome::Sys::Lock::PROXY_LOCK_PARAMS_SPEC(),
    );
    STDERR->say(sprintf 'DEBUG: Attempting to lock (in scope "%s"): %s', $self->scope, $self->resource) if $ENV{UR_DUMP_DEBUG_MESSAGES};
    my $locked = Genome::Sys::Lock->lock_resource(%params,
        resource_lock => $self->resource,
        scope => $self->scope,
    );
    if ($locked) {
        return $self;
    }
    else {
        return;
    }
}

=item unlock()

C<unlock()> is an instance method that behaves identical to
C<Genome::Sys::Lock::unlock_resource()>.

=cut

sub unlock {
    my $self = shift;
    return ! ! Genome::Sys::Lock->unlock_resource(
        resource_lock => $self->resource,
        scope => $self->scope,
    );
}

=item unlock_guard()

C<unlock_guard()> returns a Scope::Guard that will unlock the
C<Genome::Sys::LockProxy> object.

=cut

sub unlock_guard {
    my $self = shift;
    return Scope::Guard->new( sub { $self->unlock() } );
}

1;
