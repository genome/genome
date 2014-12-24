#!/usr/bin/env genome-perl

use strict;
use warnings;

BEGIN {
    $ENV{UR_DBI_NO_COMMIT} = 1;
    $ENV{UR_USE_DUMMY_AUTOGENERATED_IDS} = 1;
};

use above "Genome";
use Test::More;
use File::Compare;
use Genome::Test::Factory::SoftwareResult::User;

if (Genome::Config->arch_os ne 'x86_64') {
    plan skip_all => 'requires 64-bit machine';
}

use_ok( 'Genome::Model::Tools::DetectVariants2::Filter::PindelVafFilter');

my $refbuild_id = 101947881;
my $input_directory = $ENV{GENOME_TEST_INPUTS} . "/Genome-Model-Tools-DetectVariants2-Filter-PindelVafFilter";

my $result_users = Genome::Test::Factory::SoftwareResult::User->setup_user_hash(
    reference_sequence_build_id => $refbuild_id,
);

# Updated to v2 to allow for new columns
my $expected_dir       = $input_directory . "/expected_1/";
my $tumor_bam_file     = $input_directory . '/true_positive_tumor_validation.bam';
my $normal_bam_file    = $input_directory . '/true_positive_normal_validation.bam';
my $detector_directory = $input_directory . "/pindel-read-support-v1-";
my $test_output_dir    = Genome::Sys->create_temp_directory;

my $detector_result = Genome::Model::Tools::DetectVariants2::Result->__define__(
    output_dir       => $detector_directory,
    detector_name    => 'test',
    detector_params  => '',
    detector_version => 'awesome',
    aligned_reads    => $tumor_bam_file,
    control_aligned_reads => $normal_bam_file,
    reference_build_id    => $refbuild_id,
);
$detector_result->lookup_hash($detector_result->calculate_lookup_hash());

my $param_str;
run_test('default_params', $param_str, $result_users);

$param_str = '--variant-freq-cutoff 0.08';
run_test('non_default_params', $param_str, $result_users);

done_testing();


sub run_test {
    my ($type, $params, $result_users) = @_;
    my $output_dir = $test_output_dir."/$type";
    my $expect_dir = $expected_dir."/$type";

    my %params = (
        previous_result_id  => $detector_result->id,
        output_directory    => $output_dir,
        result_users        => $result_users,
    );

    $params{params} = $params if $params;
    my $pindel_vaf_filter = Genome::Model::Tools::DetectVariants2::Filter::PindelVafFilter->create(%params);
  
    ok($pindel_vaf_filter, "created PindelVafFilter object");
    ok($pindel_vaf_filter->execute(), "executed PindelVafFilter");

    if ($params) {
        my %parameters = split /\s+/, $params;

        for my $parameter (keys %parameters) {
            my $before_value = $parameters{$parameter};
            $parameter =~ s/^\-\-//;
            $parameter =~ s/\-/_/g;
            my $after_value  = $pindel_vaf_filter->$parameter;
            ok($before_value eq $after_value, "Parameter $parameter set correctly via params string");
        }
    }
    
    my @files = qw(indels.hq.bed indels.lq.bed);

    for my $file_name (@files) {
        my $test_output     = $output_dir."/".$file_name;
        my $expected_output = $expect_dir."/".$file_name;
        my $msg = "Output: $file_name generated as expected";
        ok(-s $test_output ,"$file_name exists and has size");
        is(compare($test_output, $expected_output), 0, $msg);
    }

    return 1;
}

