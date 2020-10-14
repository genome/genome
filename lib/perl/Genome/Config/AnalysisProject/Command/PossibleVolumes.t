#!/usr/bin/env genome-perl

use strict;
use warnings;

BEGIN {
    $ENV{UR_DBI_NO_COMMIT} = 1;
};

use Test::More tests => 2;

use above 'Genome';

use Genome::Test::Factory::AnalysisProject;

my $class = 'Genome::Config::AnalysisProject::Command::PossibleVolumes';

use_ok($class);

my $anp = Genome::Test::Factory::AnalysisProject->setup_object;

my $cmd = $class->create(analysis_project => $anp);

isa_ok($cmd, $class, 'created command');


