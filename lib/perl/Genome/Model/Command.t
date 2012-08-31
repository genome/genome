#!/usr/bin/env genome-perl

use strict;
use warnings;

use above 'Genome';

use Test::More;

use_ok('Genome::Model::Command') or die;

my $rv = eval{ Genome::Sys->shellcmd(cmd => 'genome model -h'); };
diag($@);
ok($rv, 'genome model -h');

done_testing();
exit;

