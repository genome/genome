#!/usr/bin/env genome-perl

use strict;
use warnings;

use Test::More tests => 10;

BEGIN {
    use_ok('BAP::Config');
}


my $config = BAP::Config->new();

isa_ok($config, 'BAP::Config');

my ($major, $minor, $version, $db_file ) = (
                                 $config->major(), 
                                 $config->minor(),
                                 $config->version(),
                                 $config->activity_db_file(),
                                );

ok($major =~ /^\d+$/);
ok($minor =~ /^\d+$/);
ok($major >= 2);
ok($minor >= 5);
ok($version =~ /^\d+\.\d+$/);
ok($version, join('.', $major, $minor));
ok($db_file, '/gscmnt/temp212/info/annotation/BAP_db/mgap_activity.db');

my $second_config = BAP::Config->new();

ok($config == $second_config);
