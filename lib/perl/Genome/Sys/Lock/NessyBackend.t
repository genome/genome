use strict;
use warnings;

use Genome::Sys::Lock::NessyBackend;
use Test::More;

if ($ENV{GENOME_NESSY_SERVER}) {
    plan tests => 2;
} else {
    plan skip_all => 'No Nessy URL specified for testing.';
}

use List::Util qw(shuffle);

subtest 'basic test' => sub {
    plan tests => 2;

    my $resource_name = 'NessLock.t/' . random_string();
    diag 'resource = ' . $resource_name;

    my $n = Genome::Sys::Lock::NessyBackend->new(
        url => $ENV{GENOME_NESSY_SERVER},
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
        push @r, 'NessLock.t/' . random_string();
        push @n, Genome::Sys::Lock::NessyBackend->new(
            url => $ENV{GENOME_NESSY_SERVER},
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

sub random_string {
    my @chars = map { (shuffle 'a'..'z', 'A'..'Z', 0..9)[0] } 1..10;
    return join('', @chars);
}
