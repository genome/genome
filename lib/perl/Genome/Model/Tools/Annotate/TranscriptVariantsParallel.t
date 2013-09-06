#!/usr/bin/env genome-perl

# Currently, this tests G::M::T::A::TranscriptVariantsParallel.pm by comparing its output with 
# G::M::T::A::TranscriptVariants, failing if there are any differences

use strict;
use warnings;

use Test::More skip_all => 'Test is too slow and not thorough';
use File::Temp qw/ tempdir /;
use File::Compare;
use above "Genome";

# Check test directory
my $test_dir = $ENV{GENOME_TEST_INPUTS} . '/Genome-Model-Tools-Annotate-TranscriptVariantsParallel';
my $test_output_dir = tempdir(
    '/Genome-Model-Tools-Annotate-TranscriptVariantsParallel-XXXXX',
    TMPDIR => 1,
    CLEANUP => 1);
chmod(0775,$test_output_dir);

ok (-d $test_dir, "Test data directory exists");

# Check that reference output exists in test directory
my $reference_output = $test_dir . '/transcript_variants_output.out.new';
ok (-e $reference_output, "Reference output file exists");

# Check that the test variants file exists in test directory
my $test_variants_file = $test_dir . '/variants_short.tsv';
ok (-e $test_variants_file, "Test variants file exists");

# Split by line number test
my $number_output = $test_output_dir . '/number_output';
my $number_cmd_obj = Genome::Model::Tools::Annotate::TranscriptVariantsParallel->create(
    variant_file => $test_variants_file,
    output_file => $number_output,
    split_by_number => 50,
    annotation_filter => 'top',
    cache_annotation_data_directory => 1,
);

$number_cmd_obj->execute() if $number_cmd_obj;

ok (compare($number_output, $reference_output) == 0, "Output of transcript variants and transcript variants parallel (split by line number) are the same($number_output, $reference_output)");
unlink $number_output;
#die;

# Split by chromosome test
my $chrom_output = $test_output_dir . "/chrom_output";
my $chrom_cmd_obj = Genome::Model::Tools::Annotate::TranscriptVariantsParallel->create(
    variant_file => $test_variants_file,
    output_file => $chrom_output,
    split_by_chromosome => 1,
    annotation_filter => "top",
    cache_annotation_data_directory => 1,
);

$chrom_cmd_obj->execute() if $chrom_cmd_obj;

ok (compare($chrom_output, $reference_output) == 0,  "Output of transcript variants and transcript variants parallel (split by chromosome) are the same($chrom_output, $reference_output)");
unlink $chrom_output;
