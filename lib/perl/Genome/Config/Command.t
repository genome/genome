#!/usr/bin/env genome-perl

use strict;
use warnings;

use above 'Genome';

use Test::More;

use_ok('Genome::Config::Command') or die;

my $rv = eval{ Genome::Sys->shellcmd(cmd => 'genome config -h'); };
diag($@);
ok($rv, 'genome config -h');

done_testing();
