#!/usr/bin/env genome-perl

use strict;
use warnings;
use Test::More;# 'skip_all';
use Storable 'retrieve';
use above "Genome";
use Genome::Info::UCSCConservation;

our $THIS_VERSION_ANNOTATOR_SUBCLASS = 'Genome::Model::Tools::Annotate::TranscriptVariants::Version4';

# The test variants file can hold 1..n variants
# Each variant must have a corresponding annotation
# These annotations are held in three files, one for each type of filter (top, gene, none)
# Annotations in each file are sorted according to respective filter
my ($variants, $annotations) = get_test_data();
my $ensembl_version = "67_37l_v8";
my $annotator_version = 4;
# Test annotation output for all provided variants

check_prioritization($variants, $annotations);

done_testing();


################################################################################

sub get_test_data {
    my $test_dir = $ENV{GENOME_TEST_INPUTS} . '/Genome-Transcript-VariantAnnotator/v4';
    ok (-e $test_dir, "test data directory exists at $test_dir");

    my $test_variants_file = $test_dir . "/variants.tsv";
    ok (-s $test_variants_file, "test variants file exists and has size at $test_variants_file");
    ok (-r $test_variants_file, "test variants file is readable by the user running this test $test_variants_file");

    my $none_annotations_file = $test_dir . "/none_annotations.tsv";
    ok (-s $none_annotations_file, "annotations with no filter file exists and has size");

    my $top_annotations_file = $test_dir . "/top_annotations.tsv";
    ok (-s $top_annotations_file, "annotations with top filter file exists and has size");

    my $gene_annotations_file = $test_dir . "/gene_annotations.tsv";
    ok (-s $gene_annotations_file, "annotations with gene filter file exists and has size");

    my %annotations;
    $annotations{none} = $none_annotations_file;
    $annotations{top} = $top_annotations_file;
    $annotations{gene} = $gene_annotations_file;

    return $test_variants_file, \%annotations;
}

sub check_prioritization {
    my ($variants, $annotations) = @_;
    for my $filter (qw/ none top gene /) {
        my $temp = Genome::Sys->create_temp_file_path($filter);
        my $annotator = Genome::Model::Tools::Annotate::TranscriptVariants->create(
            output_file => $temp,
            annotation_filter => $filter,
            reference_transcripts => "NCBI-human.ensembl/".$ensembl_version,
            use_version => $annotator_version,
            variant_file => $variants,
            benchmark => 1,
        );
        ok (defined $annotator, "successfully created variant annotator command");
        $annotator->dump_status_messages(1);
        ok ($annotator->execute, "successfully ran annotator command");

        compare_annotations($annotations->{$filter}, $temp);
    }
}

sub compare_annotations {
    my ($expected_annotations, $annotations) = @_;

    print "Comparing $expected_annotations to $annotations\n";

    my $diff_output = `diff $expected_annotations $annotations`;
    ok(!$diff_output, "No differences between expected output and test output") or diag("diff results:\n".$diff_output);
}
