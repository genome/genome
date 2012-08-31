#!/usr/bin/env genome-perl

use strict;
use warnings;

use Test::More tests => 4;
use File::Compare;

use above 'Genome';

BEGIN {
    use_ok('Genome::Model::Tools::Fastq::Sol2sanger');
};

my $tmp_dir = File::Temp::tempdir('Fastq-SolToSanger-XXXXX', DIR => "$ENV{GENOME_TEST_TEMP}", CLEANUP => 1);
my $fastq_file = $ENV{GENOME_TEST_INPUTS} . '/Genome-Model-Tools-Fastq/SolToSanger/test.fq';
my $expected_sanger_fastq_file = $ENV{GENOME_TEST_INPUTS} . '/Genome-Model-Tools-Fastq/SolToSanger/test.fq.sanger';

my $sol2sanger = Genome::Model::Tools::Fastq::Sol2sanger->create(
                                                                fastq_file => $fastq_file,
                                                                sanger_fastq_file => $tmp_dir .'/test.fastq',
                                                            );
isa_ok($sol2sanger,'Genome::Model::Tools::Fastq::Sol2sanger');

ok($sol2sanger->execute,'execute command '. $sol2sanger->command_name);
$DB::single = 1;
ok(compare($sol2sanger->sanger_fastq_file,$expected_sanger_fastq_file) == 0,'files are the same');

exit;
