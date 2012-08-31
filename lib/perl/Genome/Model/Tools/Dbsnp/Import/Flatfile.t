#!/usr/bin/env genome-perl

use strict;
use warnings;

BEGIN {
    $ENV{UR_DBI_NO_COMMIT} = 1;
};

use Test::More tests => 5;

use above 'Genome';

use_ok('Genome::Model::Tools::Dbsnp::Import::Flatfile');

my $sub_dir = "/Genome-Model-Tools-Dbsnp-Import-Flatfile/v2";
my $test_data_url_dir = $ENV{GENOME_TEST_URL} . $sub_dir;
my $test_data_dir = $ENV{GENOME_TEST_INPUTS} . $sub_dir;
my $test_flat_file = "$test_data_url_dir/ds_flat_chMT.txt";
my $test_output = "$test_data_dir/ds_flat_chMT.txt.out";
ok(-e $test_output, "test file $test_output exists") || die;

my $command_output = Genome::Sys->create_temp_directory();
my $output_file = $command_output."/ds_flat_chMT.bed";

my $cmd = Genome::Model::Tools::Dbsnp::Import::Flatfile->create(
    flatfile => $test_flat_file,
    output_file => $output_file,
); 

ok($cmd, 'created the importer');
ok($cmd->execute, 'importer ran successfully');
my $diff = Genome::Sys->diff_file_vs_file($test_output, $output_file);
ok(!$diff, 'returned file matches expected file')
    or diag("diff:\n" . $diff);
