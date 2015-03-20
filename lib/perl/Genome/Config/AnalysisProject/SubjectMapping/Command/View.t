#!/usr/bin/env genome-perl

use strict;
use warnings FATAL => 'all';;

use above 'Genome';
use Test::More tests => 1;

my $class = 'Genome::Config::AnalysisProject::SubjectMapping::Command::View'; 
use_ok($class);
