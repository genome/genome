#!/usr/bin/env genome-perl

use strict;
use warnings;

use above 'Genome';

use Test::More;

use_ok('Genome::Config::AnalysisProject::SubjectPairing::Command') or die;

my $rv = eval{ Genome::Sys->shellcmd(cmd => 'genome config analysis-project subject-pairing -h'); };
diag($@);
ok($rv, 'genome config analysis-project subject-pairing -h');

done_testing();
