use strict;
use warnings;

BEGIN {
    $ENV{GENOME_DEV_MODE} = 1;
};

use above 'Genome';
use Test::More tests => 3;
use Genome::Logger;

# GENOME_DEV_MODE prevent the TGI extension from loading
ok($ENV{GENOME_DEV_MODE}, 'GENOME_DEV_MODE is on');

subtest 'default logger' => sub {
    plan tests => 2;

    my $logger = Genome::Logger->logger();

    my %outputs = %{$logger->{outputs}};
    my $screen = delete $outputs{screen};

    ok($screen, 'got a screen output');
    is_deeply(\%outputs, {}, 'outputs is otherwise empty');
};

subtest 'TGI logger' => sub {
    plan tests => 3;

    require Genome::Site::TGI::Extension::Logger;
    Genome::Site::TGI::Extension::Logger->import();

    my $logger = Genome::Logger->logger();

    my %outputs = %{$logger->{outputs}};
    my $screen = delete $outputs{screen};
    my $syslog = delete $outputs{syslog};

    ok($screen, 'got a screen output');
    ok($syslog, 'got a syslog output');
    is_deeply(\%outputs, {}, 'outputs is otherwise empty');
};
