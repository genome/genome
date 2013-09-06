#!/usr/bin/env genome-perl

use strict;
use warnings;

use Test::More;
use File::Compare;
use File::Copy;
use File::Temp;

use above 'Genome';

if (`uname -a` =~ /x86_64/){
    plan tests => 4;
} else{
    plan skip_all => 'Must run on a 64 bit machine';
}
use_ok('Genome::Model::Tools::Sam::Pileup');

my $data_dir = $ENV{GENOME_TEST_INPUTS} . '/Genome-Model-Tools-Sam-Pileup';

my $tmp_dir  = File::Temp::tempdir(
    "BamToSam_XXXXXX", 
    TMPDIR => 1,
    CLEANUP => 1,
);

my $bam_file = $data_dir .'/test.bam';
my $example_output_file = $data_dir .'/expected/test.pileup.gz';

my $output_file = $tmp_dir."/test.pileup.gz";
my $refseq_path = "/gscmnt/gc4096/info/model_data/2741951221/build101947881/all_sequences.fa";
my $region_file = $data_dir."/test.region.gz";


my $pileup_cmd = Genome::Model::Tools::Sam::Pileup->create(
    bam_file => $bam_file,
    output_file => $output_file,
    samtools_version => "r599",
    reference_sequence_path => $refseq_path,
    use_bgzip => 1,
    region_file => $region_file,
);

my $result = $pileup_cmd->execute;

ok($result,"Tool exited properly.");

ok(-s $output_file,"Found output file, it has size.");

is(compare($output_file,$example_output_file), 0, "$output_file matched example pileup output file");
