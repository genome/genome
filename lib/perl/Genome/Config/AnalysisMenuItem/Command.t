#!/usr/bin/env genome-perl

use strict;
use warnings;

use above 'Genome';

use Test::More;

use_ok('Genome::Config::AnalysisMenuItem::Command') or die;

my $rv = eval{ Genome::Sys->shellcmd(cmd => 'genome config analysis-menu-item -h'); };
diag($@);
ok($rv, 'genome config analysis-menu-item -h');

done_testing();
