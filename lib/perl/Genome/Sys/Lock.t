use strict;
use warnings;

use above 'Genome';

use Genome::Sys::Lock::MockBackend;
use List::MoreUtils qw(all);
use List::Util qw(shuffle);
use Test::More tests => 10;

########################################################################
# Locking
########################################################################

subtest 'shared mandatory both lock' => sub {
    plan tests => 2;

    localize_backend_changes(sub {
        my $locks = {};
        my @backends = (
            Genome::Sys::Lock::MockBackend->new(),
            Genome::Sys::Lock::MockBackend->new(),
        );
        Genome::Sys::Lock->set_backends(site => \@backends);

        my $resource_lock = 'Lock.t/' . random_string();
        my $lock = Genome::Sys::Lock->lock_resource(
            resource_lock => $resource_lock,
            scope => 'site',
        );
        is($lock, $resource_lock, 'got the lock');
        ok((all { $_->has_lock($resource_lock) } @backends),
            'both backends locked the resource');
    });
};

subtest 'shared mandatory fails to lock' => sub {
    plan tests => 2;

    localize_backend_changes(sub {
        my $locks = {};
        my @backends = (
            Genome::Sys::Lock::MockBackend->new()
                ->locks($locks),
            Genome::Sys::Lock::MockBackend->new()
                ->locks($locks),
        );
        Genome::Sys::Lock->set_backends(site => \@backends);

        my $resource_lock = 'Lock.t/' . random_string();
        my $lock = Genome::Sys::Lock->lock_resource(
            resource_lock => $resource_lock,
            scope => 'site',
        );
        is($lock, undef, 'failed to get lock');
        ok(!(grep { $_->has_lock($resource_lock) } @backends),
            'neither backend locked the resource');
    });
};

subtest 'first optional fails to lock' => sub {
    plan tests => 3;

    localize_backend_changes(sub {
        my @backends = (
            Genome::Sys::Lock::MockBackend->new()
                ->set_always('lock', undef)
                ->set_always('is_mandatory', 0),
            Genome::Sys::Lock::MockBackend->new(),
        );

        Genome::Sys::Lock->set_backends(site => \@backends);

        my $resource_lock = 'Lock.t/' . random_string();
        my $lock = Genome::Sys::Lock->lock_resource(
            resource_lock => $resource_lock,
            scope => 'site',
        );
        is($lock, $resource_lock, 'got the lock');
        ok(!(all { $_->has_lock($resource_lock) } @backends), 'both backends did not lock the resource');
        ok((grep { $_->has_lock($resource_lock) } @backends), 'one backend did lock the resource');
    })
};

subtest 'second optional fails to lock' => sub {
    plan tests => 3;

    localize_backend_changes(sub {
        my @backends = (
            Genome::Sys::Lock::MockBackend->new(),
            Genome::Sys::Lock::MockBackend->new()
                ->set_always('lock', undef)
                ->set_always('is_mandatory', 0),
        );
        Genome::Sys::Lock->set_backends(site => \@backends);

        my $resource_lock = 'Lock.t/' . random_string();
        my $lock = Genome::Sys::Lock->lock_resource(
            resource_lock => $resource_lock,
            scope => 'site',
        );
        is($lock, $resource_lock, 'got the lock');
        ok(!(all { $_->has_lock($resource_lock) } @backends),
            'both backends did not lock the resource');
        ok((grep { $_->has_lock($resource_lock) } @backends),
            'one backend did lock the resource');
    });
};

########################################################################
# Unlocking
########################################################################

subtest 'first optional fails to lock, then unlock' => sub {
    plan tests => 3;

    localize_backend_changes(sub {
        my @backends = (
            Genome::Sys::Lock::MockBackend->new(),
            Genome::Sys::Lock::MockBackend->new()
                ->set_always('lock', undef)
                ->set_always('is_mandatory', 0),
        );
        Genome::Sys::Lock->set_backends(site => \@backends);

        my $resource_lock = 'Lock.t/' . random_string();
        my $lock = Genome::Sys::Lock->lock_resource(
            resource_lock => $resource_lock,
            scope => 'site',
        );
        is($lock, $resource_lock, 'got the lock');

        my $unlocked = Genome::Sys::Lock->unlock_resource(
            resource_lock => $resource_lock,
        );
        ok($unlocked, 'unlocked');
        ok(!(all { $_->has_lock($resource_lock) } @backends), 'both backends do not have resource lock');
    });
};

subtest 'second optional fails to lock, then unlock' => sub {
    plan tests => 3;

    localize_backend_changes(sub {
        my @backends = (
            Genome::Sys::Lock::MockBackend->new()
                ->set_always('lock', undef)
                ->set_always('is_mandatory', 0),
            Genome::Sys::Lock::MockBackend->new(),
        );
        Genome::Sys::Lock->set_backends(site => \@backends);

        my $resource_lock = 'Lock.t/' . random_string();
        my $lock = Genome::Sys::Lock->lock_resource(
            resource_lock => $resource_lock,
            scope => 'site',
        );
        is($lock, $resource_lock, 'got the lock');

        my $unlocked = Genome::Sys::Lock->unlock_resource(
            resource_lock => $resource_lock,
        );
        ok($unlocked, 'unlocked');
        ok(!(all { $_->has_lock($resource_lock) } @backends),
            'both backends do not have resource lock');
    });
};

subtest 'with both locked, first optional fails to unlock' => sub {
    plan tests => 3;

    localize_backend_changes(sub {
        my @backends = (
            Genome::Sys::Lock::MockBackend->new()
                ->set_always('unlock', undef)
                ->set_always('is_mandatory', 0),
            Genome::Sys::Lock::MockBackend->new(),
        );
        Genome::Sys::Lock->set_backends(site => \@backends);

        my $resource_lock = 'Lock.t/' . random_string();
        my $lock = Genome::Sys::Lock->lock_resource(
            resource_lock => $resource_lock,
            scope => 'site',
        );
        is($lock, $resource_lock, 'got the lock');

        my $unlocked = Genome::Sys::Lock->unlock_resource(
            resource_lock => $resource_lock,
        );
        ok($unlocked, 'unlocked');
        ok(!(all { $_->has_lock($resource_lock) } grep { $_->is_mandatory } @backends),
            'mandatory backend does not have resource lock');
    });
};

subtest 'with both locked, second optional fails to unlock' => sub {
    plan tests => 3;

    localize_backend_changes(sub {
        my @backends = (
            Genome::Sys::Lock::MockBackend->new(),
            Genome::Sys::Lock::MockBackend->new()
                ->set_always('unlock', undef)
                ->set_always('is_mandatory', 0),
        );
        Genome::Sys::Lock->set_backends(site => \@backends);

        my $resource_lock = 'Lock.t/' . random_string();
        my $lock = Genome::Sys::Lock->lock_resource(
            resource_lock => $resource_lock,
            scope => 'site',
        );
        is($lock, $resource_lock, 'got the lock');

        my $unlocked = Genome::Sys::Lock->unlock_resource(
            resource_lock => $resource_lock,
        );
        ok($unlocked, 'unlocked');
        ok(!(all { $_->has_lock($resource_lock) } grep { $_->is_mandatory } @backends),
            'mandatory backend does not have resource lock');
    });
};

subtest 'with both locked, first mandatory fails to unlock' => sub {
    plan tests => 2;

    localize_backend_changes(sub {
        my @backends = (
            Genome::Sys::Lock::MockBackend->new()
                ->set_always('unlock', undef),
            Genome::Sys::Lock::MockBackend->new(),
        );
        Genome::Sys::Lock->set_backends(site => \@backends);

        my $resource_lock = 'Lock.t/' . random_string();
        my $lock = Genome::Sys::Lock->lock_resource(
            resource_lock => $resource_lock,
            scope => 'site',
        );
        is($lock, $resource_lock, 'got the lock');

        my $unlocked = Genome::Sys::Lock->unlock_resource(
            resource_lock => $resource_lock,
        );
        ok(!$unlocked, 'unlocked failed');
    });
};

subtest 'with both locked, second mandatory fails to unlock' => sub {
    plan tests => 2;

    localize_backend_changes(sub {
        my @backends = (
            Genome::Sys::Lock::MockBackend->new(),
            Genome::Sys::Lock::MockBackend->new()
                ->set_always('unlock', undef),
        );
        Genome::Sys::Lock->set_backends(site => \@backends);

        my $resource_lock = 'Lock.t/' . random_string();
        my $lock = Genome::Sys::Lock->lock_resource(
            resource_lock => $resource_lock,
            scope => 'site',
        );
        is($lock, $resource_lock, 'got the lock');

        my $unlocked = Genome::Sys::Lock->unlock_resource(
            resource_lock => $resource_lock,
        );
        ok(!$unlocked, 'unlocked failed');
    });
};

sub random_string {
    my @chars = map { (shuffle 'a'..'z', 'A'..'Z', 0..9)[0] } 1..10;
    return join('', @chars);
}

sub localize_backend_changes {
    my $sub = shift;
    my %old_backends = Genome::Sys::Lock::all_backends();

    $sub->();

    Genome::Sys::Lock->set_backends(%old_backends);
}
