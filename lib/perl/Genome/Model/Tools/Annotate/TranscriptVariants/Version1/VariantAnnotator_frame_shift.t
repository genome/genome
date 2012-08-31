#!/usr/bin/env genome-perl

use strict;
use warnings;
use Test::More;
use above "Genome";
use File::Temp;

#This script is inteded to test the --get_frame_shift_sequence flag
my $test_dir = $ENV{GENOME_TEST_INPUTS} . '/Genome-Transcript-VariantAnnotator/get_frame_shift_sequence';
ok (-e $test_dir, "test data directory exists at $test_dir");

my $test_variants_file = $test_dir . "/variants.tsv";
ok (-s $test_variants_file, "test variants file exists and has size");

my $annotation_file = $test_dir . "/annotation_output.tsv";
ok (-s $annotation_file, "annotation file exists and has size");
my @relevant_annotation = `cat $annotation_file | grep "XM_001717859\\|NM_022552\\|NM_002520"`;
ok (scalar @relevant_annotation > 0, "successfully grabbed variants from file");

my $temp = File::Temp->new();
ok($temp, "temp file successfully created");
my $temp_filename = $temp->filename;

Genome::Model::Tools::Annotate::TranscriptVariants->execute(
    variant_file => $test_variants_file,
    reference_transcripts => "NCBI-human.combined-annotation/54_36p_v2",
    get_frame_shift_sequence => 1,
    output_file => $temp_filename, 
    annotation_filter => "none",
    use_version => 1,
);

my @relevant_new_annotation = `cat $temp_filename | grep "XM_001717859\\|NM_022552\\|NM_002520"`;
ok (scalar @relevant_new_annotation == scalar @relevant_annotation, "New annotation count matches old annotation count");

for(my $i = 1; $i < scalar @relevant_new_annotation; $i++){ #the first line is headers, so skip it
    my @new_fields = split("\t", $relevant_new_annotation[$i]);    
    my @old_fields = split("\t", $relevant_annotation[$i]); 
    my $new_aa = $new_fields[15];
    my $old_aa = $old_fields[15];
    is($new_aa, $old_aa, "$new_aa present in both annotaions");
}

done_testing();
exit;

sub variant_headers {
    return Genome::Model::Tools::Annotate->variant_attributes;
}
