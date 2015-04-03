#!/usr/bin/env genome-perl

use strict;
use warnings;

use Genome;
use Test::More tests => 3;

use Genome::Sys::Lock qw();
use Genome::Sys::Lock::MockBackend qw();
use Genome::Sys::LockMigrationProxy qw();
use Genome::Sys::LockProxy qw();
use Genome::Utility::Text qw(rand_string);
use List::Util qw(shuffle);
use Scope::Guard qw();
use Sub::Override qw();

subtest 'old fails' => sub {
    plan tests => 4;

    my %backends = (
        test_old => [ Genome::Sys::Lock::MockBackend->new() ],
        test_new => [ Genome::Sys::Lock::MockBackend->new() ],
    );
    my $backend_guard = backend_guard();
    my $scopes_guard = Sub::Override->new('Genome::Sys::Lock::scopes', sub { keys %backends });
    Genome::Sys::Lock->set_backends(%backends);

    my $resource = rand_string();
    my %lock_params = (
        block_sleep => 0,
        max_try => 0,
    );

    my $old_lock = Genome::Sys::LockProxy->new(
        resource => $resource,
        scope => 'test_old',
    )->lock(%lock_params);
    ok($old_lock, q(preemptively locked on 'test_old'));

    my $success = Genome::Sys::LockMigrationProxy->new(
        old => {
            resource => $resource,
            scope => 'test_old',
        },
        new => {
            resource => $resource,
            scope => 'test_new',
        },
    )->lock(%lock_params);
    ok(!$success, q(failed to lock when 'test_old' was already locked));

    $old_lock = Genome::Sys::LockProxy->new(
        resource => $resource,
        scope => 'test_old',
    )->lock(%lock_params);
    is($old_lock, undef, q(still locked on 'test_old'));

    my $new_lock = Genome::Sys::LockProxy->new(
        resource => $resource,
        scope => 'test_new',
    )->lock(%lock_params);
    ok($new_lock, q(able to lock on 'test_new' after LockMigrationProxy failed to lock));
};

subtest 'new fails' => sub {
    plan tests => 4;

    my %backends = (
        test_old => [ Genome::Sys::Lock::MockBackend->new() ],
        test_new => [ Genome::Sys::Lock::MockBackend->new() ],
    );
    my $backend_guard = backend_guard();
    my $scopes_guard = Sub::Override->new('Genome::Sys::Lock::scopes', sub { keys %backends });
    Genome::Sys::Lock->set_backends(%backends);

    my $resource = rand_string();
    my %lock_params = (
        block_sleep => 0,
        max_try => 0,
    );

    my $new_lock = Genome::Sys::LockProxy->new(
        resource => $resource,
        scope => 'test_new',
    )->lock(%lock_params);
    ok($new_lock, q(preemptively locked on 'test_new'));

    my $success = Genome::Sys::LockMigrationProxy->new(
        old => {
            resource => $resource,
            scope => 'test_old',
        },
        new => {
            resource => $resource,
            scope => 'test_new',
        },
    )->lock(%lock_params);
    ok(!$success, q(failed to lock when 'test_new' was already locked));

    my $old_lock = Genome::Sys::LockProxy->new(
        resource => $resource,
        scope => 'test_old',
    )->lock(%lock_params);
    ok($old_lock, q(able to lock on 'test_old' after LockMigrationProxy failed to lock));

    $new_lock = Genome::Sys::LockProxy->new(
        resource => $resource,
        scope => 'test_new',
    )->lock(%lock_params);
    is($new_lock, undef, q(still locked on 'test_new'));
};

subtest 'both ok' => sub {
    plan tests => 3;

    my %backends = (
        test_old => [ Genome::Sys::Lock::MockBackend->new() ],
        test_new => [ Genome::Sys::Lock::MockBackend->new() ],
    );
    my $backend_guard = backend_guard();
    my $scopes_guard = Sub::Override->new('Genome::Sys::Lock::scopes', sub { keys %backends });
    Genome::Sys::Lock->set_backends(%backends);

    my $resource = rand_string();
    my %lock_params = (
        block_sleep => 0,
        max_try => 0,
    );

    my $success = Genome::Sys::LockMigrationProxy->new(
        old => {
            resource => $resource,
            scope => 'test_old',
        },
        new => {
            resource => $resource,
            scope => 'test_new',
        },
    )->lock(%lock_params);
    ok($success, q(successfully locked LockMigrationProxy));

    my $old_lock = Genome::Sys::LockProxy->new(
        resource => $resource,
        scope => 'test_old',
    )->lock(%lock_params);
    is($old_lock, undef, q(failed to lock on 'test_old' after LockMigrationProxy locked));

    my $new_lock = Genome::Sys::LockProxy->new(
        resource => $resource,
        scope => 'test_new',
    )->lock(%lock_params);
    is($new_lock, undef, q(failed to lock on 'test_new' after LockMigrationProxy locked));
};

sub backend_guard {
    my %old_backends = Genome::Sys::Lock::all_backends();
    return Scope::Guard->new(sub { Genome::Sys::Lock->set_backends(%old_backends) });
}
