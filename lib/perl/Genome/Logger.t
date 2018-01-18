use strict;
use warnings;

use Test::More tests => 5;
use Test::Fatal qw(exception);

use Memoize qw();

use above 'Genome';

use_ok('Genome::Logger');

subtest 'non-color screen' => sub {
    plan tests => 1;

    no warnings 'redefine';
    local *Genome::Role::Logger::_should_color_screen = sub { 0 };
    Genome::Logger->clear_logger();

    my $logger = Genome::Logger->logger();
    isa_ok($logger->output('screen'), 'Log::Dispatch::Screen');
};

subtest 'color screen' => sub {
    plan tests => 2;

    no warnings 'redefine';
    local *Genome::Role::Logger::_should_color_screen = sub { 1 };
    Genome::Logger->clear_logger();

    ok(Genome::Logger->has_color_screen_package, 'has_color_screen_package');

    my $logger = Genome::Logger->logger();
    isa_ok($logger->output('screen'), 'Log::Dispatch::Screen::Color');
};

subtest 'context independent memoize' => sub {
    plan tests => 3;

    Genome::Logger->clear_logger();

    my $sl = Genome::Logger->logger();
    ok($sl, 'got logger in scalar context');

    my ($ll) = Genome::Logger->logger();
    ok($ll, 'got logger in list context');

    is($sl, $ll, 'both loggers are the same');
};

subtest 'exceptions' => sub {
    plan tests => 4;

    Genome::Logger->clear_logger();
    my $logger = Genome::Logger->logger();
    $logger->remove('screen');

    my $message = 'something happened';
    like(exception { Genome::Logger->croak('error', $message) },
        qr/^$message/, 'exception thrown with message');

    like(exception { Genome::Logger->croak('foo', $message) },
        qr/^invalid level/, 'exception thrown for invalid level');

    like(exception { Genome::Logger->fatal($message) },
        qr/^$message/, 'exception thrown with message');

    like(exception { Genome::Logger->fatalf($message) },
        qr/^$message/, 'exception thrown with message');
};
