#!/usr/bin/env genome-perl

use strict;
use warnings;
use Test::More;
use Storable 'retrieve';
use above "Genome";

our $THIS_VERSION_ANNOTATOR_SUBCLASS = 'Genome::Model::Tools::Annotate::TranscriptVariants::Version2';

# The test variants file can hold 1..n variants
# Each variant must have a corresponding annotation
# These annotations are held in three files, one for each type of filter (top, gene, none)
# Annotations in each file are sorted according to respective filter
my ($variants, $annotations) = get_test_data();

# Test annotation output for all provided variants
check_output($variants, $annotations->{none});

# Ensure that prioritization of annotations behaves correctly
check_prioritization($variants, $annotations);

done_testing();


################################################################################

sub variant_headers {
    return Genome::Model::Tools::Annotate->variant_attributes;
}

sub annotation_headers {
    return (
        variant_headers(),
        Genome::Model::Tools::Annotate->variant_output_attributes,
        Genome::Model::Tools::Annotate->transcript_attributes,
    );
}

sub get_test_data {
    my $test_dir = Genome::Config::get('test_inputs') . '/Genome-Transcript-VariantAnnotator';
    ok (-e $test_dir, "test data directory exists at $test_dir");

    my $test_variants_file = $test_dir . "/variants.tsv";
    ok (-s $test_variants_file, "test variants file exists and has size at $test_variants_file");
    ok (Genome::Sys->validate_file_for_reading($test_variants_file), "test variants file is readable by the user running this test $test_variants_file");

    my @variant_headers = variant_headers();
    my $variant_svr = Genome::Utility::IO::SeparatedValueReader->create(
        input => $test_variants_file,
        headers => \@variant_headers,
        separator => "\t",
        is_regex => 1,
    );
    my @variants = $variant_svr->all;
    ok (scalar @variants > 0, "successfully grabbed variants from file");

    my $none_annotations_file = $test_dir . "/none_annotations.tsv.new";
    ok (-s $none_annotations_file, "annotations with no filter file exists and has size");

    my $top_annotations_file = $test_dir . "/top_annotations.tsv.new";
    ok (-s $top_annotations_file, "annotations with top filter file exists and has size");

    my $gene_annotations_file = $test_dir . "/gene_annotations.tsv.new";
    ok (-s $gene_annotations_file, "annotations with gene filter file exists and has size");

    my %annotations;
    $annotations{none} = retrieve($none_annotations_file);
    $annotations{top} = retrieve($top_annotations_file);
    $annotations{gene} = retrieve($gene_annotations_file);

    ok (scalar @{$annotations{none}} > 0, "succesfully grabbed test annotations for filter \'none\'");
    ok (scalar @{$annotations{top}} > 0, "successfully grabbed test annotations for filter \'top\'");
    ok (scalar @{$annotations{gene}} > 0, "successfully grabbed test annotations for filter \'gene\'");

    return \@variants, \%annotations;
}

sub create_annotator {
    my $annotation_model = Genome::Model->get(name => 'NCBI-human.combined-annotation');
    my $annotation_build = $annotation_model->build_by_version('54_36p_v2');

    my @data_directories = $annotation_build->determine_data_directory();
    my $data_directory;
    if (@data_directories < 2) {
        $data_directory = $data_directories[0];
    } else {
        $data_directory = \@data_directories;
    }

    my $annotator = $THIS_VERSION_ANNOTATOR_SUBCLASS->create(
        data_directory => $data_directory,
    );
    return $annotator;
}

sub get_annotations_for_variant {
    my ($variant, $annotations) = @_;
    my @annotations_for_variant;

    for my $anno (@$annotations) {
        if ($anno->{chromosome_name} eq $variant->{chromosome_name} and
            $anno->{start} eq $variant->{start} and
            $anno->{stop} eq $variant->{stop}) 
        {
            push @annotations_for_variant, $anno;
        }
    }

    return @annotations_for_variant;
}

sub check_output {
    my ($variants, $annotations) = @_;
    my $variant_num = 0;

    my $output_annotator = create_annotator();
    ok (defined $output_annotator, "succesfully created variant annotator object");

    for my $variant (@$variants) {
        $variant->{type} = Genome::Model::Tools::Annotate->infer_variant_type($variant);
        my @expected_annotations = get_annotations_for_variant($variant, $annotations);

        my @annotations = $output_annotator->transcripts(%$variant);
        # The expected data also contains the keys/values from the variant hash
        # Go ahead and insert the variant data into each annotation
        foreach (@annotations) {
            %$_ = ( %$_, %$variant);
        }

        ok (compare_annotations(\@expected_annotations, \@annotations), "annotation output matches expected output for variant $variant_num");
        $variant_num++;
    }
}

sub check_prioritization {
    my ($variants, $annotations) = @_;
    for my $filter (qw/ none top gene /) {
        my $annotator = create_annotator();
        my $variant_num = 0;
        for my $variant (@$variants) {
            my @test_output = get_annotations_for_variant($variant, $annotations->{$filter});
            my @output;
            if ($filter eq "none") {
                @output = $annotator->transcripts(%$variant);
            }
            elsif ($filter eq "gene") {
                @output = $annotator->prioritized_transcripts(%$variant);
            }
            elsif ($filter eq "top") {
                my $output = $annotator->prioritized_transcript(%$variant);
                push @output, $output;
            }

            for my $out (@output) {
                for my $key (keys %$variant) {
                    $out->{$key} = $variant->{$key};
                }
            }

            is(scalar(@output), scalar(@test_output), "received same number of results as expected for variant $variant_num with filter $filter");
            ok (compare_annotations(\@output, \@test_output), "annotation ordering matches expected after prioritization for variant $variant_num with filter $filter");
            $variant_num++;
        }
    }
}

sub compare_annotations {
    my ($expected_annotations, $annotations) = @_;

    my @expected_annotations = sort { $a->{'transcript_name'} cmp $b->{'transcript_name'} }
                               @$expected_annotations;

    my @annotations = sort { $a->{'transcript_name'} cmp $b->{'transcript_name'} }
                      @$annotations;

    for (my $i = 0; $i < scalar @expected_annotations; $i++) {
        my $expected = $expected_annotations[$i];
        my $annotation = $annotations[$i];
        for my $key (sort keys %$expected) {
            unless ($expected->{$key} eq $annotation->{$key}) {
                diag "*** for transcript name ". $expected->{'transcript_name'}
                     . " attribute $key : expected-> " . $expected->{$key} . " received-> " . $annotation->{$key} . " ***\n";
                return 0;
            }
        }
    }
    return 1;
}
