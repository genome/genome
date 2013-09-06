#!/usr/bin/env genome-perl

use strict;
use warnings;

use Test::More tests => 7;
use File::Compare;

use above 'Genome';

BEGIN {
    use_ok('Genome::Model::Tools::Fastq::Trimq2::Fragment');
};

my $base_dir = $ENV{GENOME_TEST_INPUTS} . '/Genome-Model-Tools-Fastq/Trimq2/Fragment';

my $tmp_dir = File::Temp::tempdir(
    'Fastq-Trimq2-XXXXX', 
    TMPDIR => 1,
    CLEANUP => 1
);
my $fastq_file = "$base_dir/test_fragment.fastq";

#Test sanger fastq, report yes, length_limit 32 bp
my $trimq2 = Genome::Model::Tools::Fastq::Trimq2::Fragment->create(
    fastq_file  => $fastq_file,
    output_dir  => $tmp_dir,
    #trim_string => '#',  default is #
);
isa_ok($trimq2,'Genome::Model::Tools::Fastq::Trimq2');

ok($trimq2->execute,'execute command '. $trimq2->command_name);

is($trimq2->out_file, "$tmp_dir/test_fragment.trimq2.fastq", 'output name is ok');

for my $file qw(test_fragment.trimq2.fastq trimq2.fragment.filtered.fastq trimq2.report) {
    my $output_file = $tmp_dir."/$file";
    my $expect_file = $base_dir."/$file";
    ok(compare($output_file, $expect_file) == 0, "Output $file is created as expected");
}

#TODO should also add Illumina fastq trim test
