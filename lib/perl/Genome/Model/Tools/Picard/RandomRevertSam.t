#!/usr/bin/env genome-perl

use strict;
use warnings;

use above 'Genome';

use Test::More tests => 4;
use File::Compare;

use_ok('Genome::Model::Tools::Picard::RandomRevertSam');

my $data_dir = $ENV{GENOME_TEST_INPUTS} . '/Genome-Model-Tools-Picard-RandomRevertSam';
my $test_bam = $data_dir .'/test.bam';
my $expected_bam = $data_dir .'/expected.bam';

my $tmp_dir = File::Temp::tempdir('Picard-RandomRevertSam-'.Genome::Sys->username.'-XXXX',DIR => "$ENV{GENOME_TEST_TEMP}",CLEANUP => 1);
my $tmp_bam = $tmp_dir .'/tmp.bam';

my $rand_revert_sam = Genome::Model::Tools::Picard::RandomRevertSam->create(
    input_file => $test_bam,
    output_file => $tmp_bam,
    probability => '0.5',
    random_seed => '123',
    remove_duplicate_information => 1,
    remove_alignment_information => 1,
);
isa_ok($rand_revert_sam,'Genome::Model::Tools::Picard::RandomRevertSam');
ok($rand_revert_sam->execute,'execute command '. $rand_revert_sam->command_name);
ok(!compare($expected_bam,$tmp_bam),'Expected BAM '. $expected_bam .' is identical to the RandomRevertSam BAM file '. $tmp_bam);


exit;
