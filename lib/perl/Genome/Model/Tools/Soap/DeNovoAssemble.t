#!/usr/bin/env genome-perl

use strict;
use warnings;

use above "Genome";
use Test::More;

my $archos = `uname -a`;
unless ($archos =~ /64/) {
    plan skip_all => "Must run from 64-bit machine";
}

#check for test data files
my $data_dir = $ENV{GENOME_TEST_INPUTS} . '/Genome-Model-Tools-Soap/DeNovoAssemble/data_dir';
ok(-d $data_dir, "Data dir exists");

my @data_files = qw/ config.txt
                     SRR038746.pair1.fastq SRR038746.pair2.fastq SRR038746.single.fastq
                     SRR042027.pair1.fastq SRR042027.pair2.fastq SRR042027.single.fastq /;

foreach (@data_files) {
    ok(-s $data_dir."/$_", "Data dir $_ file exists");
}

#make temp test dir
my $temp_dir = Genome::Sys->create_temp_directory();
ok(-d $temp_dir, "Temp test dir created");

#run soap denovo
my $create = Genome::Model::Tools::Soap::DeNovoAssemble->create(
    version => 1.04,
    config_file => "$data_dir/config.txt",
    kmer_size => 31,
    resolve_repeats => 1,
    kmer_frequency_cutoff => 1,
    cpus => 8,
    output_dir_and_file_prefix => "$temp_dir/TEST",
    );

ok( $create, "Created gmt soap de-novo assemble");
ok(($create->execute) == 1, "Command ran successfully");

#compare output files
my $data_run_dir = $ENV{GENOME_TEST_INPUTS} . '/Genome-Model-Tools-Soap/DeNovoAssemble/run_dir';
ok(-d $data_run_dir, "Data run dir exists");

#soap-denovo output files
my @run_files = qw/ TEST.Arc          TEST.contig         TEST.gapSeq        TEST.links       TEST.newContigIndex
                    TEST.peGrads      TEST.preGraphBasic  TEST.readOnContig  TEST.scafSeq     TEST.updated.edge
                    TEST.ContigIndex  TEST.edge           TEST.kmerFreq      TEST.markOnEdge  TEST.path
                    TEST.preArc       TEST.readInGap      TEST.scaf          TEST.scaf_gap    TEST.vertex /;

foreach (@run_files) {
    ok(-s $data_run_dir."/$_", "Data run dir $_ exists");
}

my @diffs = `diff -r --brief $temp_dir $data_run_dir`;
is (scalar (@diffs), 0, "Run outputs match");

#<STDIN>;

done_testing();

exit;
