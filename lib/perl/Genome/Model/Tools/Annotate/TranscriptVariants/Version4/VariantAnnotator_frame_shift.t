#!/usr/bin/env genome-perl

use strict;
use warnings;
use Test::More;
use above "Genome";
use File::Temp;

#This script is inteded to test the --get_frame_shift_sequence flag
my $test_dir = $ENV{GENOME_TEST_INPUTS} . '/Genome-Transcript-VariantAnnotator/get_frame_shift_sequence/v4';
ok (-e $test_dir, "test data directory exists at $test_dir");

my $test_variants_file = $test_dir . "/variants.tsv";
ok (-s $test_variants_file, "test variants file exists and has size");

my $annotation_file = $test_dir . "/annotation_output.tsv";
ok (-s $annotation_file, "annotation file exists and has size");
my @relevant_annotation = `cat $annotation_file | grep "ENST00000301067\\|ENST00000342783\\|ENST00000375547\\|ENST00000563486"`;
ok (scalar @relevant_annotation > 0, "successfully grabbed variants from file");

my $temp = Genome::Sys->create_temp_file_path;
ok($temp, "temp file successfully created");

Genome::Model::Tools::Annotate::TranscriptVariants->execute(
    variant_file => $test_variants_file,
    reference_transcripts => "NCBI-human.ensembl/67_37l_v8",
    get_frame_shift_sequence => 1,
    output_file => $temp, 
    annotation_filter => "none",
    use_version => 4,
);

my @relevant_new_annotation = `cat $temp | grep "ENST00000301067\\|ENST00000342783\\|ENST00000375547\\|ENST00000563486"`;
ok (scalar @relevant_new_annotation == scalar @relevant_annotation, "New annotation count matches old annotation count");

for(my $i = 1; $i < scalar @relevant_new_annotation; $i++){ #the first line is headers, so skip it
    my @new_fields = split("\t", $relevant_new_annotation[$i]);    
    my @old_fields = split("\t", $relevant_annotation[$i]); 
    my $new_aa = $new_fields[15];
    my $old_aa = $old_fields[15];
    is($new_aa, $old_aa, "$new_aa present in both annotaions");
}

done_testing();
