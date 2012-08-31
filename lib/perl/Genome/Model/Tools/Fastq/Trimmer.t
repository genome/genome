#!/usr/bin/env genome-perl

use strict;
use warnings;

use above 'Genome';

use Test::More tests => 4;
use File::Compare;

use_ok('Genome::Model::Tools::Fastq::Trimmer');

my $tmp_dir = File::Temp::tempdir('Fastq-Trimmer-XXXXX', DIR => "$ENV{GENOME_TEST_TEMP}", CLEANUP => 1);

my $data_dir = $ENV{GENOME_TEST_INPUTS} . '/Genome-Model-Tools-Fastq-Trimmer';
#my $tmp_dir = $data_dir;
my $expected_fastq = $data_dir .'/5_25.fastq';
my $fastq_file = $data_dir .'/original.fastq';


my $trim_5_25 = Genome::Model::Tools::Fastq::Trimmer->create(
    input_fastq => $fastq_file,
    output_fastq => $tmp_dir .'/tmp-5_25.fastq',
    five_prime => 5,
    three_prime => 25
);
isa_ok($trim_5_25,'Genome::Model::Tools::Fastq::Trimmer');
ok($trim_5_25->execute,'execute command '. $trim_5_25->command_name);
ok(!compare($trim_5_25->output_fastq,$expected_fastq),'expected fastq');

exit;
