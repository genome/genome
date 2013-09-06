#!/usr/bin/env genome-perl

use strict;
use warnings;

use above "Genome";
use Test::More tests => 5;
use File::Compare;
use File::Temp qw(tempfile);
use File::Path qw(rmtree);

use_ok( 'Genome::Model::Tools::Breakdancer::TigraToBreakdancer');

my $test_input_dir  = $ENV{GENOME_TEST_INPUTS} . '/Genome-Model-Tools-Breakdancer-TigraToBreakdancer/';
my $sv_file     = $test_input_dir . 'breakdancer.out';
my $tigra_file  = $test_input_dir . 'tigra.out';
my $pass_file   = $test_input_dir . 'pass.out';
my $fail_file   = $test_input_dir . 'fail.out';

my $tmp_dir = File::Temp::tempdir(
    'Genome-Model-Tools-Breakdancer-TigraToBreakdancer-XXXXX', 
    TMPDIR => 1,
    CLEANUP => 1,
);

my($tmp_pass_file, $tmp_fail_file) = map{$tmp_dir.'/'.$_}qw(tmp_pass.out tmp_fail.out);

my $tigra_2_bd = Genome::Model::Tools::Breakdancer::TigraToBreakdancer->create(
    original_breakdancer_file => $sv_file,
    tigra_output_file         => $tigra_file,
    pass_filter_file          => $tmp_pass_file,
    fail_filter_file          => $tmp_fail_file,
);

ok($tigra_2_bd, 'TigraToBreakdancer created ok');
ok($tigra_2_bd->execute(), 'TigraToBreakdancer executed ok');

is(compare($pass_file, $tmp_pass_file), 0, 'pass_filter_file output as expected');
is(compare($fail_file, $tmp_fail_file), 0, 'fail_filter_file output as expected');
