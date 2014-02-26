#!/usr/bin/env genome-perl
use strict;
use warnings;
use above "Genome";
use Test::More;
use Genome::Command::Tester qw/run_and_diff/;

run_and_diff(results_version => '2014-02-25');
done_testing();

# This infers the module to test from the test name,
# and infers the command to run from the synopsis
# and replaces /tmp/output_dir with an autogeerated temp directory.
# To change results, set a new version/date, and run with "REBUILD".

