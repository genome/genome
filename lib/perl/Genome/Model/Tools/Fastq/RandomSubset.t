#!/usr/bin/env genome-perl

use strict;
use warnings;

use above 'Genome';

use Test::More tests => 11;
use File::Compare;

use_ok('Genome::Model::Tools::Fastq::RandomSubset');

my $tmp_dir = File::Temp::tempdir('Fastq-RandomSubset-XXXXX', DIR => "$ENV{GENOME_TEST_TEMP}", CLEANUP => 1);

my $data_dir = $ENV{GENOME_TEST_INPUTS} . '/Genome-Model-Tools-Fastq/RandomSubset';
#my $tmp_dir = $data_dir;
my $expected_50_file = $data_dir .'/50_seq.fastq';
my $expected_10_file = $data_dir .'/10_seq.fastq';
my $expected_100bp_file = $data_dir .'/100bp_seq.fastq';
my $fastq_file = $data_dir .'/All.fastq';


my $rs_50 = Genome::Model::Tools::Fastq::RandomSubset->create(
        input_read_1_fastq_files => [$fastq_file],
        output_read_1_fastq_file => $tmp_dir .'/tmp-50.fastq',
        #output_fastq_file => $expected_50_file,
        limit_type => 'reads',
        limit_value => 50,
        seed_phrase => 'test_seed',
        );
isa_ok($rs_50,'Genome::Model::Tools::Fastq::RandomSubset');
ok($rs_50->execute,'execute command '. $rs_50->command_name);
ok(!compare($rs_50->output_read_1_fastq_file,$expected_50_file),'expected 50 file equal');

my $rs_10 = Genome::Model::Tools::Fastq::RandomSubset->create(
        input_read_1_fastq_files => [$fastq_file],
        output_read_1_fastq_file => $tmp_dir .'/tmp-10.fastq',
        #output_fastq_file => $expected_10_file,
        limit_type => 'reads',
        limit_value => 10,
        seed_phrase => 'test_seed',
        );
isa_ok($rs_10,'Genome::Model::Tools::Fastq::RandomSubset');
ok($rs_10->execute,'execute command '. $rs_10->command_name);
ok(!compare($rs_10->output_read_1_fastq_file,$expected_10_file),'expected 10 file equal');
is($rs_10->_shortest_seq,36,'36 base pair sequence length');

my $rs_100bp = Genome::Model::Tools::Fastq::RandomSubset->create(
        input_read_1_fastq_files => [$fastq_file],
        output_read_1_fastq_file => $tmp_dir .'/tmp-100bp.fastq',
        #output_fastq_file => $expected_100bp_file,
        limit_type => 'base_pair',
        limit_value => 100,
        seed_phrase => 'test_seed',
        );
isa_ok($rs_100bp,'Genome::Model::Tools::Fastq::RandomSubset');
ok($rs_100bp->execute,'execute command '. $rs_100bp->command_name);
ok(!compare($rs_100bp->output_read_1_fastq_file,$expected_100bp_file),'expected 100bp file equal');
exit;
