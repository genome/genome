#!/usr/bin/env genome-perl

use strict;
use warnings FATAL => 'all';;

use above 'Genome';
use Test::More tests => 1;

use_ok('Genome::Config::AnalysisProject::Command::View');
