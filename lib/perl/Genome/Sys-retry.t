use strict;
use warnings;

use above 'Genome';
use Test::More tests => 4;
use Test::Exception;

subtest 'retry when return undef' => sub {
    plan tests => 2;

    my $tries = 3;
    my $executions = 0;
    my $callback = sub {
        $executions++;
        return;
    };
    my $rv = Genome::Sys::retry(callback => $callback, tries => $tries, delay => 0);
    is($executions, $tries, sprintf(q(executed %d times), $tries));
    ok(!$rv, 'returned falsey');
};

subtest 'retry when returning 0' => sub {
    plan tests => 2;

    my $tries = 3;
    my $executions = 0;
    my $callback = sub {
        $executions++;
        return 0;
    };
    my $rv = Genome::Sys::retry(callback => $callback, tries => $tries, delay => 0);
    is($executions, $tries, sprintf(q(executed %d times), $tries));
    ok(!$rv, 'returned falsey');
};

subtest 'retry when returning ($executions %% $stop_at == 0)' => sub {
    plan tests => 2;

    my $tries = 3;
    my $executions = 0;
    my $stop_at = 2;
    my $callback = sub {
        $executions++;
        return ($executions % $stop_at == 0);
    };
    my $rv = Genome::Sys::retry(callback => $callback, tries => $tries, delay => 0);
    is($executions, $stop_at, sprintf(q(executed %d times), $stop_at));
    ok($rv, 'returned truthy');
};

subtest 'type validation' => sub {
    plan tests => 7;

    throws_ok { Genome::Sys::retry(tries => 1, delay => 0) }
        qr/'callback'/, 'absent callback throws exception';

    throws_ok { Genome::Sys::retry(callback => sub {}, delay => 0) }
        qr/'tries'/, 'absent tries throws exception';

    throws_ok { Genome::Sys::retry(callback => sub {}, tries => 1) }
        qr/'delay'/, 'absent delay throws exception';

    throws_ok { Genome::Sys::retry(callback => sub {}, tries => 0, delay => 0) }
        qr/is greater than zero/, 'tries => 0 throws exception';

    throws_ok { Genome::Sys::retry(callback => sub {}, tries => 1.5, delay => 0) }
        qr/is an integer/, 'tries => 1.5 throws exception';

    throws_ok { Genome::Sys::retry(callback => sub {}, tries => 0, delay => -1) }
        qr/is greater than zero/, 'delay => -1 throws exception';

    throws_ok { Genome::Sys::retry(callback => sub {}, tries => 1, delay => 1.5) }
        qr/is an integer/, 'delays => 1.5 throws exception';
};
