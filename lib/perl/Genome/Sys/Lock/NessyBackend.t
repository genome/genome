use strict;
use warnings;

use Genome::Sys::Lock::NessyBackend;
use Genome::Utility::Text qw(rand_string);
use Test::More;

my $nessy_server = Genome::Config::get('nessy_server');
if ($nessy_server) {
    plan tests => 3;
} else {
    plan skip_all => 'No Nessy URL specified for testing.';
}

use List::Util qw(shuffle);

subtest 'basic test' => sub {
    plan tests => 2;

    my $resource_name = 'NessLock.t/' . rand_string();
    diag 'resource = ' . $resource_name;

    my $n = Genome::Sys::Lock::NessyBackend->new(
        url => $nessy_server,
        is_mandatory => 1,
    );

    my $resource = $n->lock(
        resource => $resource_name,
        timeout => 5,
        wait_announce_interval => 10,
    );
    is($resource, $resource_name, 'locked a Nessy resource');

    my $unlocked = $n->unlock($resource);
    is($unlocked, 1, 'unlocked a Nessy resource');
};

subtest 'instance validation' => sub {
    plan tests => 4;

    my (@r, @n);
    for (1..2) {
        push @r, 'NessLock.t/' . rand_string();
        push @n, Genome::Sys::Lock::NessyBackend->new(
            url => $nessy_server,
            is_mandatory => 1,
        );
        $n[-1]->lock(
            resource => $r[-1],
            timeout => 5,
            wait_announce_interval => 10,
        );
    }

    isnt($n[0]->client, $n[1]->client, 'instances do not have the same client');
    ok(!$n[0]->has_lock($r[1]), 'first instance does not have lock on second resource');
    ok(!$n[1]->has_lock($r[0]), 'second instance does not have lock on first resource');

    my $lock = $n[1]->lock(
        resource => $r[0],
        timeout => 5,
        wait_announce_interval => 10,
    );
    ok(!$lock, 'second client failed to get lock on first resource (and did not crash)');

    $n[0]->unlock($r[0]);
    $n[1]->unlock($r[1]);
};

subtest 'namespace validation' => sub {
    plan tests => 7;

    my $resource = rand_string();
    my $root_n = Genome::Sys::Lock::NessyBackend->new(
        url => $nessy_server,
        is_mandatory => 1,
    );
    my $namespace_n = Genome::Sys::Lock::NessyBackend->new(
        url => $nessy_server,
        is_mandatory => 1,
        namespace => 'foo',
    );

    isnt($root_n->client, $namespace_n->client, 'instances do not have the same client');
    ok(!$root_n->has_lock($resource), 'root instance does not have lock on resource');
    ok(!$namespace_n->has_lock($resource), 'namespace instance does not have lock on resource');

    my $root_lock = $root_n->lock(
        resource => $resource,
        timeout => 5,
        wait_announce_interval => 10,
    );
    ok($root_lock, 'root instance got lock on resource');

    my $namespace_lock = $namespace_n->lock(
        resource => $resource,
        timeout => 5,
        wait_announce_interval => 10,
    );
    ok($namespace_lock, 'namespace instance got lock on resource');

    ok($root_n->unlock($resource), 'root instance unlocked');
    ok($namespace_n->unlock($resource), 'namespace instance unlocked');
};
