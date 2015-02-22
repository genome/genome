#!/usr/bin/env genome-perl

use strict;
use warnings;

use above 'Genome';
use Test::More tests => 31;

use_ok('Genome::Model::Tools::Vcf::VcfFilter');

my $test_dir = $ENV{GENOME_TEST_INPUTS} . "/Genome-Model-Tools-Vcf-VcfFilter";

# V2 adds the FT per-sample column
# V3 adds the FT header line mistakenly omitted from v2
# V4 corrects the FT header to be a string
# V5 makes this test actually meaningful. The previous test data was broken. It set bed_input but snvs.hq was not a bed file.
# V6 changes the expected output file name to snvs and indels since we are now testing indels
# V7 change the description of format "FT" to match that of TCGA-compliant format
my $expected_base = "expected.v7";
my $expected_dir = "$test_dir/$expected_base";

for my $input_type ("unfiltered_input", "filtered_input") {
    for my $variant_type ("snvs", "indels") {
        my $output_file = Genome::Sys->create_temp_file_path;
        my $input_dir = "$test_dir/$input_type";
        my $input_vcf = "$input_dir/$variant_type.vcf.gz";
        my $hq_filter_file = "$input_dir/$variant_type.hq.bed";
        my $expected_file = "$expected_dir/$variant_type.vcf";

        my $command= Genome::Model::Tools::Vcf::VcfFilter->create(
            output_file => $output_file,
            vcf_file => $input_vcf,
            filter_file => $hq_filter_file,
            filter_keep => 1,
            filter_name => "TEST",
            filter_description => "TEST",
            bed_input => 1,
        );

        ok($command, 'Command created');
        my $rv = $command->execute;
        ok($rv, 'Command completed successfully');
        ok(-s $output_file, "output file created");

        # The files will have a timestamp that will differ. Ignore this but check the rest.
        my $diff = Genome::Sys->diff_file_vs_file($output_file, $expected_file);
        ok(!$diff, 'output matched expected result')
            or diag("diff results:\n" . $diff);

        # FIXME is this crazy? This would have caught bad test data the first time. But if it is manually checked when tests are updated it would be fine...
        my $wc = `wc -l $hq_filter_file`;
        my ($expected_pass_count) = split /\s/, $wc;
        my $actual_pass_count = `grep :PASS $output_file | wc -l`;
        chomp($actual_pass_count);

        # FIXME this is pretty hacky for a test... but two lines in the VCF do not match the bed file because they are complex (RPL) deletions.
        if ($variant_type eq "indels") {
            $expected_pass_count = $expected_pass_count - 2;
        }

        is($actual_pass_count, $expected_pass_count, "Number of passing variants is as expected");
    }
}

for my $input_type ("large_input", "somatic_with_null_sample_input") {
    my $variant_type = "snvs";
    my $output_file = Genome::Sys->create_temp_file_path;
    my $input_dir = "$test_dir/$input_type";
    my $input_vcf = "$input_dir/$variant_type.vcf.gz";
    my $hq_filter_file = "$input_dir/$variant_type.hq.bed";
    my $expected_file = "$input_dir/output.vcf";

    my $command= Genome::Model::Tools::Vcf::VcfFilter->create(
        output_file => $output_file,
        vcf_file => $input_vcf,
        filter_file => $hq_filter_file,
        filter_keep => 1,
        filter_name => "TEST",
        filter_description => "TEST",
        bed_input => 1,
    );

    ok($command, 'Command created');
    my $rv = $command->execute;
    ok($rv, 'Command completed successfully');
    ok(-s $output_file, "output file created");

    # The files will have a timestamp that will differ. Ignore this but check the rest.
    my $diff = Genome::Sys->diff_file_vs_file($output_file, $expected_file);
    system "diff $output_file $expected_file";
    ok(!$diff, 'output matched expected result')
        or diag("diff results:\n" . $diff);

    # FIXME is this crazy? This would have caught bad test data the first time. But if it is manually checked when tests are updated it would be fine...
    my $wc = `wc -l $hq_filter_file`;
    my ($expected_pass_count) = split /\s/, $wc;
    my $actual_pass_count = `grep :PASS $output_file | wc -l`;
    chomp($actual_pass_count);

    # FIXME this is pretty hacky for a test... but two lines in the VCF do not match the bed file because they are complex (RPL) deletions.
    if ($variant_type eq "indels") {
        $expected_pass_count = $expected_pass_count - 2;
    }

    is($actual_pass_count, $expected_pass_count, "Number of passing variants is as expected");
}
