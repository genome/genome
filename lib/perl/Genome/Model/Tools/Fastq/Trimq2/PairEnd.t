#!/usr/bin/env genome-perl

use strict;
use warnings;

use Test::More tests => 12;
use File::Compare;

use above 'Genome';

BEGIN {
    use_ok('Genome::Model::Tools::Fastq::Trimq2::PairEnd');
};

my $base_dir = $ENV{GENOME_TEST_INPUTS} . '/Genome-Model-Tools-Fastq/Trimq2/PairEnd';

my $tmp_dir = File::Temp::tempdir(
    'Fastq-Trimq2-XXXXX', 
    TMPDIR => 1,
    CLEANUP => 1
);

my $pe1_fq = $base_dir.'/test_pair_end_1.fastq';
my $pe2_fq = $base_dir.'/test_pair_end_2.fastq';

#Test sanger fastq, report yes, length_limit 32 bp
my $sanger_trimq2 = Genome::Model::Tools::Fastq::Trimq2::PairEnd->create(
    pair1_fastq_file => $pe1_fq,
    pair2_fastq_file => $pe2_fq,
    output_dir       => $tmp_dir,
    trim_string      => 'E',
);
isa_ok($sanger_trimq2,'Genome::Model::Tools::Fastq::Trimq2');

ok($sanger_trimq2->execute,'execute command '. $sanger_trimq2->command_name);

is($sanger_trimq2->pair1_out_file, "$tmp_dir/test_pair_end_1.trimq2.fastq", 'pair_end 1 output name ok');
is($sanger_trimq2->pair2_out_file, "$tmp_dir/test_pair_end_2.trimq2.fastq", 'pair_end 2 output name ok');
is($sanger_trimq2->pair_as_frag_file, "$tmp_dir/trimq2.pair_as_fragment.fastq", 'pair as fragment output name ok');

for my $file qw(trimq2.pair_as_fragment.fastq
test_pair_end_1.trimq2.fastq test_pair_end_2.trimq2.fastq trimq2.pair_end.filtered.fastq trimq2.pair_as_fragment.filtered.fastq trimq2.report) {
    my $output_file = $tmp_dir."/$file";
    my $expect_file = $base_dir."/$file";
    ok(compare($output_file, $expect_file) == 0, "Output $file is created as expected");
}

#TODO should also add Illumina fastq trim test
