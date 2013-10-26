#!/usr/bin/env genome-perl

use strict;
use warnings;

use above 'Genome';
use File::Compare;
use Test::More;

# This test was auto-generated because './Model/Tools/SamStat/HtmlReport.pm'
# had no '.t' file beside it.  Please remove this test if you believe it was
# created unnecessarily.  This is a bare minimum test that just compiles Perl
# and the UR class.
use_ok('Genome::Model::Tools::SamStat::HtmlReport');

my $dir = $ENV{GENOME_TEST_INPUTS} . '/Genome-Model-Tools-SamStat-HtmlReport/v1';
my $tmp_dir = Genome::Sys->create_temp_directory(Genome::Sys->username . "SamStat_XXXXXX");

my $bam_link = $dir . '/test.bam';
my $bam_file = readlink $bam_link;

my $test_link = $tmp_dir . '/test.bam';
Genome::Sys->create_symlink($bam_file, $test_link);

my $test_out     = $test_link . '.html';
my $expected_out = $bam_link . '.html';

my $cmd = Genome::Model::Tools::SamStat::HtmlReport->create(
    input_files => $test_link,
    fix_xp      => 0,
);

ok($cmd->execute, 'samstat htmlreport executed ok');
cmp_ok(compare($test_out, $expected_out), '==', 0, "samstat html created ok.");

done_testing();
