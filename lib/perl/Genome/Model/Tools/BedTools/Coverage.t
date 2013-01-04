#!/usr/bin/env genome-perl

use strict;
use warnings;

use Test::More skip_all => 'broken pending review by jwalker'; 
#plan tests => 5;
use File::Compare;

use above 'Genome';

use_ok('Genome::Model::Tools::BedTools');
use_ok('Genome::Model::Tools::BedTools::Coverage');

my $tmp_dir = File::Temp::tempdir('BedTools-Coverage-'.Genome::Sys->username.'-XXXX',DIR => "$ENV{GENOME_TEST_TEMP}",CLEANUP => 1);

my $data_dir = $ENV{GENOME_TEST_INPUTS} . '/Genome-Model-Tools-BedTools-Coverage';

my $bam_file = $data_dir .'/test.bam';
my $regions_file = $data_dir .'/test_regions.bed';
my $expected_output_file = $data_dir .'/test.out';

my $cov = Genome::Model::Tools::BedTools::Coverage->create(
    output_file => $tmp_dir .'/test.out',
    input_file_a => $bam_file,
    input_file_b => $regions_file,
    histogram => 1,
);
isa_ok($cov,'Genome::Model::Tools::BedTools::Coverage');
ok($cov->execute,'execute Coverage command '. $cov->command_name);

ok(!compare($expected_output_file,$cov->output_file),'expected coverage file '. $expected_output_file .' is identical to '. $cov->output_file);
