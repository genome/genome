#!/usr/bin/env genome-perl

use strict;
use warnings;

use above "Genome";

use Test::More tests => 1;

use_ok('Genome::Qc::Command::Config::Diff') or die;

done_testing;
