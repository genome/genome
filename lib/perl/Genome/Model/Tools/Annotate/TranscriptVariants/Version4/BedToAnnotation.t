#!/usr/bin/env genome-perl

use strict;
use warnings;
use Test::More;
use File::Temp qw(tempfile);
use above "Genome";

our $THIS_VERSION_ADAPTOR_SUBCLASS = 'Genome::Model::Tools::Annotate::TranscriptVariants::Version4::BedToAnnotation';

my $test_dir = $ENV{GENOME_TEST_INPUTS} . '/Genome-Model-Tools-Annotate-TranscriptVariants-BedToAnnotation';
ok (-e $test_dir, "test data directory exists at $test_dir");

test_snvs($test_dir);
test_indels($test_dir);
test_extra_columns($test_dir);

done_testing();

sub test_snvs{
    my $test_dir = shift;

    my $test_bed_file = $test_dir . "/snvs.bed";
    ok (-s $test_bed_file, "snv bed file exists and has size at $test_bed_file");
    ok (-r $test_bed_file, "snv bed file is readable by the user running this test $test_bed_file");

    my $annotation_file = $test_dir. "/snvs.annotation";
    ok (-s $annotation_file, "snv output file exists and has size at $annotation_file");
    ok (-r $annotation_file, "snv variants file is readable by the user running this test $annotation_file");

    my ($output_fh, $output_file) = tempfile(UNLINK => 1);
    my $adaptor = $THIS_VERSION_ADAPTOR_SUBCLASS->create(
        snv_file => $test_bed_file,
        output => $output_file,
    );

    ok(defined $adaptor, "created adaptor for snv test");
    ok($adaptor->execute, "executed snv adaptor");
    my $diff_output = `diff $output_file $annotation_file`;
    ok(!$diff_output, "snv adaptor output diffed clean");
    $output_fh->close;
}


sub test_indels{
    my $test_dir = shift;

    my $test_bed_file = $test_dir . "/indels.bed";
    ok (-s $test_bed_file, "indel bed file exists and has size at $test_bed_file");
    ok (-r $test_bed_file, "indel bed file is readable by the user running this test $test_bed_file");

    my $annotation_file = $test_dir. "/indels.annotation";
    ok (-s $annotation_file, "indel output file exists and has size at $annotation_file");
    ok (-r $annotation_file, "indel variants file is readable by the user running this test $annotation_file");

    my ($output_fh, $output_file) = tempfile(UNLINK => 1);
    my $adaptor = $THIS_VERSION_ADAPTOR_SUBCLASS->create(
        indel_file => $test_bed_file,
        output => $output_file,
    );

    ok(defined $adaptor, "created adaptor for indel test");
    ok($adaptor->execute, "executed indel adaptor");
    my $diff_output = `diff $output_file $annotation_file`;
    print "$diff_output\n";
    ok(!$diff_output, "indel adaptor output diffed clean");
    $output_fh->close;
}


sub test_extra_columns{
    my $test_dir = shift;

    my $test_bed_file = $test_dir . "/extra_columns.bed";
    ok (-s $test_bed_file, "extra_column bed file exists and has size at $test_bed_file");
    ok (-r $test_bed_file, "extra_column bed file is readable by the user running this test $test_bed_file");

    my $annotation_file = $test_dir. "/extra_columns.annotation";
    ok (-s $annotation_file, "extra_column output file exists and has size at $annotation_file");
    ok (-r $annotation_file, "extra_column variants file is readable by the user running this test $annotation_file");

    my ($output_fh, $output_file) = tempfile(UNLINK => 1);
    my $adaptor = $THIS_VERSION_ADAPTOR_SUBCLASS->create(
        snv_file => $test_bed_file,
        output => $output_file,
        extra_columns => "type,score",
    );

    ok(defined $adaptor, "created adaptor for extra_column test");
    ok($adaptor->execute, "executed extra_column adaptor");
    my $diff_output = `diff $output_file $annotation_file`;
    ok(!$diff_output, "extra_column adaptor output diffed clean");
    $output_fh->close;
}
