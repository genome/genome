#!/usr/bin/env genome-perl

use strict;
use warnings;

use Test::More tests => 4;
use File::Compare;

use above 'Genome';

BEGIN {
    use_ok('Genome::Model::Tools::Fastq::Sol2phred');
};

# prints solexa ascii characters for 0-40
#perl -e 'for (0 .. 40) { print  chr($_ + 64) ."\n"; }'
# prints phred ascii characters for 0-40
#perl -e 'for (0 .. 40) { print  chr($_ + 33) ."\n"; }'

#For the new Solexa-Phred qualities with an offset of 64, the equation
#simplifies to
#  $fastq = chr(ord($solq) - 64 + 33);
#or just
#  $fastq = chr(ord($solq) - 31);

my $tmp_dir = File::Temp::tempdir('Fastq-SolToPhred-XXXXX', DIR => "$ENV{GENOME_TEST_TEMP}", CLEANUP => 1);
my $fastq_file = $ENV{GENOME_TEST_INPUTS} . '/Genome-Model-Tools-Fastq/SolToPhred/test.fq';
my $expected_phred_fastq_file = $ENV{GENOME_TEST_INPUTS} . '/Genome-Model-Tools-Fastq/SolToPhred/test.fq.utf8.phred';

my $sol2phred = Genome::Model::Tools::Fastq::Sol2phred->create(
                                                                fastq_file => $fastq_file,
                                                                phred_fastq_file => $tmp_dir .'/test.fastq',
                                                            );
isa_ok($sol2phred,'Genome::Model::Tools::Fastq::Sol2phred');

ok($sol2phred->execute,'execute command '. $sol2phred->command_name);
ok(compare($sol2phred->phred_fastq_file,$expected_phred_fastq_file) == 0,'files are the same');
