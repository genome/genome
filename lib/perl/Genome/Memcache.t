#!/usr/bin/env genome-perl


use above "Genome";
use Data::Dumper;

use Test::More tests => 4;


use_ok('Genome::Memcache');

SKIP: {
    skip "should rewrite test to not depend on production memcached server", 3;

    my $m = Genome::Memcache->server();
    ok($m , 'server object');

    my $key = 'memcache.t';
        if (my $then = $m->get($key)) {
            diag('deleteing previous test key from: ' . $then);
            $m->delete($key);
        }

    my $now = UR::Context->current->now();
    ok($m->set($key, $now, 10), 'setting cache item');

    cmp_ok($m->get($key), 'eq', $now, "getting cache item $now");
}

done_testing;
