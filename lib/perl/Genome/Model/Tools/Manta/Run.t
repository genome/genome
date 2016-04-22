#!/usr/bin/env genome-perl

use strict;
use warnings;

use above 'Genome';
use File::Spec;
use Test::More tests => 2;

my $pkg = 'Genome::Model::Tools::Manta::Run';
use_ok ($pkg);

my $tmp_dir = Genome::Sys->create_temp_directory();
my $cpus = 4;
my @expected_cmdline = (
   File::Spec->join($tmp_dir,$pkg->_tool_subcommand_name),
   '--jobs',
   $cpus,
   '--mode',
   'local',
);

my $run_manta = Genome::Model::Tools::Manta::Run->create(
   working_directory => $tmp_dir,
   jobs => $cpus,
);

my $cmdline_array_ref = $run_manta->build_cmdline_array_ref();
is_deeply($cmdline_array_ref,\@expected_cmdline,'expected command-line array ref');
