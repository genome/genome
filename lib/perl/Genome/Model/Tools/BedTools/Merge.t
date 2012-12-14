#!/usr/bin/env genome-perl

use strict;
use warnings;

use Test::More tests => 5;
use File::Compare;

use above 'Genome';

use_ok('Genome::Model::Tools::BedTools');
use_ok('Genome::Model::Tools::BedTools::Merge');

my $tmp_dir = File::Temp::tempdir('BedTools-Merge-'.Genome::Sys->username.'-XXXX',DIR => "$ENV{GENOME_TEST_TEMP}",CLEANUP => 1);

my $data_dir = $ENV{GENOME_TEST_INPUTS} . '/Genome-Model-Tools-BedTools-Merge';

my $bed_file = $data_dir .'/full.bed';
my $expected_file = $data_dir .'/merged.bed';

my $merge = Genome::Model::Tools::BedTools::Merge->create(
    output_file => $tmp_dir .'/test.bed',
    input_file => $bed_file,
);
isa_ok($merge,'Genome::Model::Tools::BedTools::Merge');
ok($merge->execute,'execute mergeBed command '. $merge->command_name);

ok(!compare($expected_file,$merge->output_file),'expected BED file '. $expected_file .' is identical to '. $merge->output_file);
