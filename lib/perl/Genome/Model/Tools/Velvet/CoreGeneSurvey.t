#!/usr/bin/env genome-perl

use strict;
use warnings;

use Test::More;
use above "Genome";
require File::Compare;

use_ok ('Genome::Model::Tools::Velvet::CoreGeneSurvey') or die;
use_ok ('Genome::Model::Tools::Bacterial::CoreGeneCoverage') or die;

done_testing();
