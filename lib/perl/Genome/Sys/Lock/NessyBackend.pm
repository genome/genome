package Genome::Sys::Lock::NessyBackend;

use strict;
use warnings;

use Carp qw(carp croak);
use POSIX qw(strftime);
use Sys::Hostname qw(hostname);

use Genome::Logger;

use Mouse;
with qw(Genome::Sys::Lock::Backend);

has 'url' => (is => 'ro', isa => 'Str');
has 'client' => (is => 'ro', isa => 'Nessy::Client', lazy_build => 1);
has 'claims' => (is => 'rw', isa => 'ArrayRef', auto_deref => 1);

# TTL is specified in seconds.
my $DEFAULT_TTL = 2 * 60 * 60;


sub lock {
    my ($self, %args) = @_;

    my $resource = $args{resource};
    unless (defined $resource) {
        croak('resource not defined');
    }

    my $timeout = delete $args{timeout};
    unless (defined $timeout) {
        croak('timeout not defined');
    }

    my $wait_announce_interval = delete $args{wait_announce_interval};
    unless (defined $wait_announce_interval) {
        croak('wait_announce_interval not defined');
    }

    return unless $self->client;

    if ($self->_is_holding_nessy_lock($resource)) {
        Genome::Logger->fatal("Tried to lock resource more than once: $resource");
    }

    my $claim_warning = '';
    my $initial_time = time();

    my %user_data = (
        host => hostname,
        pid => $$,
        lsf_id => ($ENV{LSB_JOBID} || 'NONE'),
        user => scalar(getpwuid($<)),
        genome_build_id => ($ENV{GENOME_BUILD_ID} || 'NONE'),
        lsf_project => ($ENV{WF_LSF_PROJECT} || 'NONE'),
        requested_at => strftime('%a, %d %b %Y %T %z', localtime($initial_time)),
    );
    my $info_content = join("\n", map { $_ . ': ' . $user_data{$_} } keys %user_data);

    my $wait_announce_timer = AnyEvent->timer(
        after => $wait_announce_interval,
        interval => $wait_announce_interval,
        cb => sub {
            my $total_elapsed_time = time() - $initial_time;
            Genome::Logger->notice("waiting (total_elapsed_time = $total_elapsed_time seconds) on lock for resource '$resource': $claim_warning. lock_info is:\n$info_content");
        },
    );
    my $claim = $self->client->claim($resource, timeout => $timeout,
        ttl => $DEFAULT_TTL, user_data => \%user_data);
    undef $wait_announce_timer;
    if ($claim) {
        $self->add_claim($resource => $claim);
        return $resource;
    }
    return;
}

sub unlock {
    my ($self, $resource) = @_;
    unless ($resource) {
        carp('resource is not set');
    }

    if ($self->has_client) {
        my $claim = $self->claim($resource);
        $self->remove_claim($claim);
        if ($claim) {
            $claim->release;
        } else {
            Genome::Logger->error("Nessy tried to release, but no claim in slot for resource: $resource");
        }
    } else {
        return 1;
    }
}

sub release_all {
    my $self = shift;

    for my $claim ( $self->claims ) {
        my $resource = $claim->resource_name;
        warn("Removing remaining lock: '$resource'") unless $ENV{'HARNESS_ACTIVE'};
        $self->unlock($resource);
    }
}

sub translate_lock_args {
    my ($self, %args) = @_;

    my $block_sleep = delete $args{block_sleep};

    $args{timeout} = $self->_new_style_lock_timeout_from_args(
        block_sleep => $block_sleep,
        max_try     => delete $args{max_try},
    );

    $args{wait_announce_interval} = _translate_wait_announce_interval(
        block_sleep => $block_sleep,
        wait_announce_interval => delete $args{wait_announce_interval},
    );

    $args{resource} = delete $args{resource_lock};

    return %args;
}

sub translate_unlock_args {
    my ($self, %args) = @_;
    return delete $args{resource_lock};
}

sub is_enabled {
    my $self = shift;
    return $self->url ? 1 : 0;
}

sub _build_client {
    my $self = shift;

    if ($self->url) {
        require Nessy::Client;
        return Nessy::Client->new(url => $self->url);
    }

    return;
}

sub add_claim {
    my ($self, $resource, $claim) = @_;
    $self->claims([$self->claims, $claim]);
}

sub remove_claim {
    my ($self, $claim) = @_;
    my @claims = grep { $_ ne $claim } $self->claims;
    $self->claims(\@claims);
}

sub clear_claims {
    my ($self) = @_;
    $self->claims([]);
}

sub claim {
    my ($self, $resource) = @_;
    my ($claim) = grep { $_->resource_name eq $resource } $self->claims;
    return $claim;
}

sub has_lock { _is_holding_nessy_lock(@_) }
sub _is_holding_nessy_lock {
    my ($self, $resource) = @_;
    return $self->claim($resource) ? 1 : 0;
}

sub min_timeout {
    return 5;
}

sub _new_style_lock_timeout_from_args {
    my ($self, %args) = @_;

    my $block_sleep = delete $args{block_sleep} || 0;

    my $max_try = delete $args{max_try} || 0;

    my $min_timeout = min_timeout();
    my $timeout = $max_try * $block_sleep;
    unless ($timeout >= $min_timeout) {
        $timeout = $min_timeout;
        carp("increasing timeout to minimum ($min_timeout)");
    }

    return $timeout;
}

sub _translate_wait_announce_interval {
    my (%args) = @_;

    my $block_sleep = delete $args{block_sleep};
    my $wait_announce_interval = delete $args{wait_announce_interval};

    if ($wait_announce_interval) {
        return $wait_announce_interval;
    }

    if ($block_sleep) {
        return $block_sleep;
    }

    croak 'cannot translate wait_announce_interval';
}

__PACKAGE__->meta->make_immutable();
