use strict;
use warnings;

use above 'Genome';
use Test::More tests => 3;

subtest 'retry when return undef' => sub {
    plan tests => 2;

    my $retries = 3;
    my $executions = 0;
    my $callback = sub {
        $executions++;
        return;
    };
    my $rv = Genome::Sys::retry(callback => $callback, retries => $retries, delay => 0);
    is($executions, $retries, sprintf(q(executed %d times), $retries));
    ok(!$rv, 'returned falsey');
};

subtest 'retry when returning 0' => sub {
    plan tests => 2;

    my $retries = 3;
    my $executions = 0;
    my $callback = sub {
        $executions++;
        return 0;
    };
    my $rv = Genome::Sys::retry(callback => $callback, retries => $retries, delay => 0);
    is($executions, $retries, sprintf(q(executed %d times), $retries));
    ok(!$rv, 'returned falsey');
};

subtest 'retry when returning ($executions %% $stop_at == 0)' => sub {
    plan tests => 2;

    my $retries = 3;
    my $executions = 0;
    my $stop_at = 2;
    my $callback = sub {
        $executions++;
        return ($executions % $stop_at == 0);
    };
    my $rv = Genome::Sys::retry(callback => $callback, retries => $retries, delay => 0);
    is($executions, $stop_at, sprintf(q(executed %d times), $stop_at));
    ok($rv, 'returned truthy');
};
