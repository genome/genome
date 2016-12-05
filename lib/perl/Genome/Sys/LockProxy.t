#!/usr/bin/env genome-perl

use strict;
use warnings;

use above 'Genome';
use Test::More tests => 5;

use Genome::Sys::Lock qw();
use Genome::Sys::LockProxy qw();
use Genome::Sys::Lock::MockBackend qw();
use Genome::Utility::Text qw(rand_string);
use List::Util qw(shuffle);
use Scope::Guard qw();
use Sub::Override qw();

my %backends = (
    site => [ Genome::Sys::Lock::MockBackend->new() ],
);
my $backend_guard = backend_guard();
Genome::Sys::Lock->set_backends(%backends);

my $lock = Genome::Sys::LockProxy->new(
    resource => rand_string(),
    scope => 'site',
);

ok($lock->lock, 'lock succeeds');
ok(!$lock->lock, 'second lock fails');
ok($lock->unlock, 'unlock succeeds');
ok($lock->lock, 'lock after unlock succeeds');
ok($lock->unlock, 'final unlock succeeds');

sub backend_guard {
    my %old_backends = Genome::Sys::Lock::all_backends();
    return Scope::Guard->new(sub { Genome::Sys::Lock->set_backends(%old_backends) });
}
