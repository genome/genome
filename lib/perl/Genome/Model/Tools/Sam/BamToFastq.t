#!/usr/bin/env genome-perl

use strict;
use warnings;

use Test::More;
use File::Compare;

use above 'Genome';

if (`uname -a` =~ /x86_64/){
    plan tests => 19;
} else{
    plan skip_all => 'Must run on a 64 bit machine';
}

use_ok('Genome::Model::Tools::Sam::BamToFastq');

my $tmp_dir = Genome::Sys->create_temp_directory('Genome-Model-Tools-Sam-BamToFastq-'. Genome::Sys->username);
my $data_dir = $ENV{GENOME_TEST_INPUTS} . '/Genome-Model-Tools-Sam-BamToFastq';
my $bam_file = $data_dir .'/test.bam';

my $i = 0;
#All Reads
&run_tests(expected_reads => 1000);

#All Unmapped Reads
&run_tests(
    include_flag => '0x0008',
    expected_reads => 65,
);

#All Mapped Reads
&run_tests(
    exclude_flag => '0x0008',
    expected_reads => 935,
);

#All Duplicate Reads
&run_tests(
    include_flag => '0x0400',
    expected_reads => 346,
);

#All Properly Paired
&run_tests(
    include_flag => '0x0002',
    expected_reads => 753,
);

#All Unproperly Paired
&run_tests(
    exclude_flag => '0x0002',
    expected_reads => 247,
);

#TODO: Test combining flags to get say all mapped duplicates read 1

exit;

sub run_tests {
    my %params = @_;
    my $expected_reads = delete($params{expected_reads});
    my $fastq_file = $tmp_dir .'/'. $i++ .'.fastq';
    $params{fastq_file} = $fastq_file;
    $params{bam_file} = $bam_file;
    my $cmd = Genome::Model::Tools::Sam::BamToFastq->create(%params);
    isa_ok($cmd,'Genome::Model::Tools::Sam::BamToFastq');
    ok($cmd->execute,'execute BamToFastq command '. $cmd->command_name);
    my $reads = count_reads_in_fastq_file($fastq_file);
    is($reads,$expected_reads,$expected_reads .' reads converted from BAM to FASTQ');
}

sub count_reads_in_fastq_file {
    my $fastq_file = shift;
    my $wc = `wc -l $fastq_file`;
    unless ($wc) {
        die('Unable to count reads in FASTQ files (No output from `wc` returned.): '. $fastq_file);
    }
    my ($lines) = $wc =~ /^(\d+)/;
    if ($lines % 4) {
        die('FASTQ file '. $fastq_file .' has a line count of '. $lines .' which is not divisible by four!');
    }
    my $reads = $lines/4; #four lines per read in FASTQ file format
    return $reads;
}
