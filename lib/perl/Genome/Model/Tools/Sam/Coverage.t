#!/usr/bin/env genome-perl

use strict;
use warnings;

use above "Genome";

use Test::More;
use File::Temp;
use File::Copy;
use File::Compare;

if (Genome::Config->arch_os ne 'x86_64') {
    plan skip_all => 'requires 64-bit machine';
}
else {
    plan tests => 4;
}

use_ok('Genome::Model::Tools::Sam::Coverage');

my $data_dir = $ENV{GENOME_TEST_INPUTS} . '/Genome-Model-Tools-Sam-Coverage';

my $tmp_dir  = File::Temp::tempdir(
    "Coverage_XXXXXX", 
    TMPDIR => 1,
    CLEANUP => 1,
);

my $compare_to_file = $ENV{GENOME_TEST_INPUTS} . '/Genome-Model-Tools-Sam-Coverage/compare.txt';

my $aligned_file_name = "normal.tiny.bam"; 
my $output_file_name = "coverage.out";

my $ref_file = Genome::Config::reference_sequence_directory() . '/NCBI-human-build36/all_sequences.fa';

my $aligned_file = $data_dir."/".$aligned_file_name;
my $output_file = $tmp_dir."/".$output_file_name;

print "Input aligned reads file: $aligned_file\n";
print "Ref file: $ref_file\n";
print "Output file: $output_file\n";

my $coverage = Genome::Model::Tools::Sam::Coverage->create(
    aligned_reads_file => $aligned_file,
    reference_file => $ref_file,
    output_file => $output_file,                                                      
    return_output => 1,
    #use_version => 'r350wu1',  #This only works from r350wu1, once set Sam.pm default_value to r350wu1, this will become unnecessary
    coverage_command => $ENV{GENOME_SW} . '/samtools/bamcheck/bamcheck-v0.13/bam-check -q 1',
);

isa_ok($coverage,'Genome::Model::Tools::Sam::Coverage');
my $result = $coverage->execute;
$result =~ m/\nAverage Coverage:(\S+)/g; 
my $haploid_coverage=$1 if defined($1);
ok( $haploid_coverage eq '0.000', "haploid coverage calculated correctly" );

cmp_ok(compare($output_file, $compare_to_file), '==', 0, 'Coverage output file create ok.');

#This test runs but has been commented out to save time in the test suite.
#my $coverage2 = Genome::Model::Tools::Sam::Coverage->create(
#    aligned_reads_file => $aligned_file,
#    reference_file => $ref_file,
#    return_output => 1,                                                     
#);

#my $result = $coverage2->execute;
#print "Result: \n".$result;
#ok(defined($result),'return value defined');

#old result string
#$result =~ m/Average depth across all non-gap regions: (\S+)/g; 

#$result =~ m/Average Coverage:(\S+)/g; 
#my $haploid_coverage=$1 if defined($1);

#ok( $haploid_coverage eq '223.681', "haploid coverage calculated correctly" );
