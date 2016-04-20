#!/usr/bin/env genome-perl

use strict;
use warnings;

use above 'Genome';
use Test::More tests => 2;

my $pkg = 'Genome::Model::Tools::Manta::Run';
use_ok ($pkg);

my $expected_cmdline = '/tmp/runWorkflow.py --mode "local"';

my $run_manta = Genome::Model::Tools::Manta::Run->create(
   working_directory => '/tmp',
);

my $cmdline_string = $run_manta->build_cmdline_string();
is($cmdline_string,$expected_cmdline,'expected command-line');

exit;