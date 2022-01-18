#!/usr/bin/env genome-perl

use strict;
use warnings;
use Test::More;
use File::Temp qw(tempfile);
use above "Genome";

our $THIS_VERSION_ADAPTOR_SUBCLASS = 'Genome::Model::Tools::Annotate::TranscriptVariants::Version2::BedToAnnotation';

my $test_dir = Genome::Config::get('test_inputs') . '/Genome-Model-Tools-Annotate-TranscriptVariants-BedToAnnotation';
ok (-e $test_dir, "test data directory exists at $test_dir");

test_snvs($test_dir);
test_indels($test_dir);
test_nonsense($test_dir);

done_testing();

sub test_snvs{
    my $test_dir = shift;

    my $test_bed_file = $test_dir . "/snvs.bed";
    ok (-s $test_bed_file, "snv bed file exists and has size at $test_bed_file");
    ok (Genome::Sys->validate_file_for_reading($test_bed_file), "snv bed file is readable by the user running this test $test_bed_file");

    my $annotation_file = $test_dir. "/snvs.annotation";
    ok (-s $annotation_file, "snv output file exists and has size at $annotation_file");
    ok (Genome::Sys->validate_file_for_reading($annotation_file), "snv variants file is readable by the user running this test $annotation_file");

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
    ok (Genome::Sys->validate_file_for_reading($test_bed_file), "indel bed file is readable by the user running this test $test_bed_file");

    my $annotation_file = $test_dir. "/indels.annotation";
    ok (-s $annotation_file, "indel output file exists and has size at $annotation_file");
    ok (Genome::Sys->validate_file_for_reading($annotation_file), "indel variants file is readable by the user running this test $annotation_file");

    my ($output_fh, $output_file) = tempfile(UNLINK => 1);
    my $adaptor = $THIS_VERSION_ADAPTOR_SUBCLASS->create(
        indel_file => $test_bed_file,
        output => $output_file,
    );

    ok(defined $adaptor, "created adaptor for indel test");
    ok($adaptor->execute, "executed indel adaptor");
    my $diff_output = `diff $output_file $annotation_file`;
    ok(!$diff_output, "indel adaptor output diffed clean");
    $output_fh->close;
}

sub test_nonsense {
    my $test_dir = shift;

    my $test_bed_file = $test_dir . '/nonsense.bed';
    ok(-s $test_bed_file, "nonsense bed file exists and has size at $test_bed_file");
    ok(Genome::Sys->validate_file_for_reading($test_bed_file), "nonsense bed file is readable by user running this test $test_bed_file");

    my $annotation_file = $test_dir. "/nonsense.annotation";
    ok (-z $annotation_file, "nonsense output file exists and has 0 size at $annotation_file");
    ok (Genome::Sys->validate_file_for_reading($annotation_file), "nonsense variants file is readable by the user running this test $annotation_file");

    my ($output_fh, $output_file) = tempfile(UNLINK => 1);
    my $adaptor = $THIS_VERSION_ADAPTOR_SUBCLASS->create(
        indel_file => $test_bed_file,
        output => $output_file,
    );

    ok(defined $adaptor, "created adaptor for nonsense test");
    ok($adaptor->execute, "executed nonsense adaptor");
    my $diff_output = `diff $output_file $annotation_file`;
    ok(!$diff_output, "nonsense adaptor output diffed clean");
    $output_fh->close;
}
