use strict;
use warnings;

use Test::More tests => 4;

use Memoize qw();

use_ok('Genome::Logger');

subtest 'non-color screen' => sub {
    plan tests => 1;

    no warnings 'redefine';
    local *Genome::Logger::should_color_screen = sub { 0 };
    Memoize::flush_cache('Genome::Logger::logger');

    my $logger = Genome::Logger->logger();
    isa_ok($logger->output('screen'), 'Log::Dispatch::Screen');
};

subtest 'color screen' => sub {
    plan tests => 2;

    no warnings 'redefine';
    local *Genome::Logger::should_color_screen = sub { 1 };
    Memoize::flush_cache('Genome::Logger::logger');

    ok(Genome::Logger->has_color_screen_package, 'has_color_screen_package');

    my $logger = Genome::Logger->logger();
    isa_ok($logger->output('screen'), 'Log::Dispatch::Screen::Color');
};

subtest 'context independent memoize' => sub {
    plan tests => 3;

    my $sl = Genome::Logger->logger();
    ok($sl, 'got logger in scalar context');

    my ($ll) = Genome::Logger->logger();
    ok($ll, 'got logger in list context');

    is($sl, $ll, 'both loggers are the same');
};
