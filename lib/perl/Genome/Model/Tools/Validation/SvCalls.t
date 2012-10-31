#!/usr/bin/env genome-perl

use strict;
use warnings;
use above 'Genome';
use File::Compare;
use Test::More;
use Cwd;

if (Genome::Config->arch_os ne 'x86_64') {
    plan skip_all => 'requires 64-bit machine';
}
else {
    plan tests => 31;
}

#test that module is use-able
use_ok('Genome::Model::Tools::Validation::SvCalls');

#set up input filenames
my $test_dir = $ENV{GENOME_TEST_INPUTS} . "/Genome-Model-Tools-Validation-SvCalls/v2";
my @input_filenames = qw(bd_del_call sd_del_call test.t.9.bam test.n.9.bam test.t.9.bam.bai test.n.9.bam.bai);
@input_filenames = map { "$test_dir/$_" } @input_filenames;

#set up temp dir to run tests, and copy input files over (not optimal, but the modules used write full paths of input files to the headers)
my $temp_dir = Genome::Sys->create_temp_directory() . "/";
for my $input (@input_filenames) { system("cp $input $temp_dir"); }

my $input_bd_call = $test_dir . "bd_del_call";
my $input_sd_call = $test_dir . "sd_del_call";
my $input_t_bam = $test_dir . "test.t.9.bam";
my $input_n_bam = $test_dir . "test.n.9.bam";

#create an array for the expected output files
my @output_filenames = qw(test.assembly_input test.assembly_output.csv test.assembly_output.cm test.assembly_output.fasta test.assembly_output.csv.index test.assembly_output.fasta.merged test.assembly_output.csv.merged test.assembly_output.csv.merged.readcounts test.assembly_output.csv.merged.readcounts.noevents test.assembly_output.csv.merged.readcounts.germline test.assembly_output.csv.merged.readcounts.somatic test.assembly_output.csv.merged.readcounts.ambiguous test.assembly_output.csv.merged.readcounts.somatic.wgs_readcounts test.assembly_output.csv.merged.readcounts.somatic.wgs_readcounts.somatic);
my @expected_output_files = map { "$test_dir/$_" } @output_filenames;

#create an array for the current test output files
my @test_output_files = map { "$temp_dir/$_" } @output_filenames;

#set up test command in temp dir
my $wd = cwd(); chdir($temp_dir);
my $sv_calls = Genome::Model::Tools::Validation::SvCalls->create(
    assembled_call_files => 'bd_del_call',
    squaredancer_files => 'sd_del_call',
    normal_val_bam => 'test.n.9.bam',
    normal_wgs_bam => 'test.n.9.bam',
    tumor_val_bam => 'test.t.9.bam',
    tumor_wgs_bam => 'test.t.9.bam',
    output_filename_prefix => "test",
    patient_id => 'TEST',
    build => '36',
);

#test that object was created successfully
ok($sv_calls, 'Created SvCalls object');

#test that object was executed successfully
ok($sv_calls->execute(), 'Executed SvCalls object OK');

#run file comparisons
for my $i (0..$#test_output_files) {
    ok(-s $test_output_files[$i], 'Generated output file: ' . $output_filenames[$i] . ' ok');
    is(compare($test_output_files[$i], $expected_output_files[$i]), 0, 'Output matched expected results for ' . $output_filenames[$i]);
}

#change back to original working dir
chdir($wd);
