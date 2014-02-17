#!/usr/bin/env genome-perl

use strict;
use warnings;

use above 'Genome';

use Test::More;

use_ok('Genome::Config::AnalysisProject::SubjectMapping::Command') or die;

my $rv = eval{ Genome::Sys->shellcmd(cmd => 'genome config analysis-project subject-mapping -h'); };
diag($@);
ok($rv, 'genome config analysis-project subject-mapping -h');

done_testing();
