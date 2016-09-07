#!/usr/bin/env genome-perl

use strict;
use warnings;

use Test::More;

if (Genome::Sys->arch_os ne 'x86_64') {
    plan skip_all => 'requires 64-bit machine';
}

use above 'Genome';
use Genome::File::Vcf::Genotype;
use Genome::Utility::Test qw(compare_ok);

use_ok('Genome::Model::Tools::Bed::Convert::VcfToBed');

my $instance = Genome::Model::Tools::Bed::Convert::VcfToBed->create( source => "fake", output => "faker");
ok(defined $instance, "able to create an instance of the object");

subtest 'SNV gt het conversion' => sub {
    my $ref = 'G';
    my @var = ('G','T');
    my $result = $instance->_convert_snv_gt_to_bed($ref, @var);
    is_deeply($result, ['G', 'K'], "heterozygous conversion as expected");
};

subtest 'SNV gt hom conversion' => sub {
    my $ref = 'G';
    my @var = ('T','T');
    my $result = $instance->_convert_snv_gt_to_bed($ref, @var);
    is_deeply($result, ['G', 'T'], "homozygous conversion as expected");
};

subtest 'SNV gt biallelic conversion' => sub {
    my $ref = 'G';
    my @var = ('T','C');

    my $result = $instance->_convert_snv_gt_to_bed($ref, @var);
    is_deeply($result, ['G', 'Y'], "bi-allelic het conversion as expected");
};

my $base_test_dir = Genome::Config::get('test_inputs') . "/Genome-Model-Tools-Bed-Convert-Vcf-To-Bed/v2/";
my $temp_base_dir = File::Temp::tempdir('VcfToBedXXXXX', CLEANUP => 1, TMPDIR => 1);

subtest 'Indel Integration test' => sub {
    my $indel_vcf_input = $base_test_dir . "test_indel.vcf";
    my $expected_bed = $base_test_dir . "expected_indel.bed";
    my $test_output = $temp_base_dir . "/" . "test_indel.bed";

    my $converter = Genome::Model::Tools::Bed::Convert::VcfToBed->create(
       source => $indel_vcf_input,
       output => $test_output,
    );

    ok($converter, 'converter created');
    my $rv = $converter->execute;
    is($rv, 1, 'Testing for succesful execution. Expecting 1. Got : ' . $rv);

    ok(-s $test_output, "output file created");

    compare_ok($expected_bed, $test_output, name => 'output matched expected results');
};

subtest 'Indel Integration test one based' => sub {
    my $indel_vcf_input = $base_test_dir . "test_indel.vcf";
    my $expected_bed = $base_test_dir . "expected_indel_one_based.bed";
    my $test_output = $temp_base_dir . "/" . "test_indel_one_based.bed";

    my $converter = Genome::Model::Tools::Bed::Convert::VcfToBed->create(
       source => $indel_vcf_input,
       output => $test_output,
       one_based => 1,
    );

    ok($converter, 'converter created');
    my $rv = $converter->execute;
    is($rv, 1, 'Testing for succesful execution. Expecting 1. Got : ' . $rv);

    ok(-s $test_output, "output file created");

    compare_ok($expected_bed, $test_output, name => 'output matched expected results');
};

subtest 'SNV Integration test' => sub {
    my $snv_vcf_input = $base_test_dir . "test_snv.vcf";
    my $expected_bed = $base_test_dir . "expected_snv.bed";
    my $test_output = $temp_base_dir . "/" . "test_snv.bed";

    my $converter = Genome::Model::Tools::Bed::Convert::VcfToBed->create(
       source => $snv_vcf_input,
       output => $test_output,
    );

    ok($converter, 'converter created');
    my $rv = $converter->execute;
    is($rv, 1, 'Testing for succesful execution. Expecting 1. Got : ' . $rv);

    ok(-s $test_output, "output file created");

    compare_ok($expected_bed, $test_output, name => 'output matched expected results');

};

subtest 'SNV multisample integration test' => sub {
    my $snv_vcf_input = $base_test_dir . "test_snv.vcf";
    my $expected_bed = $base_test_dir . "expected_snv2.bed";
    my $test_output = $temp_base_dir . "/" . "test_snv2.bed";

    my $converter = Genome::Model::Tools::Bed::Convert::VcfToBed->create(
       source => $snv_vcf_input,
       sample_name => "H_NJ-HCC1395-HCC1395_BL",
       output => $test_output,
    );

    ok($converter, 'converter created');
    my $rv = $converter->execute;
    is($rv, 1, 'Testing for succesful execution. Expecting 1. Got : ' . $rv);

    ok(-s $test_output, "output file created");

    compare_ok($expected_bed, $test_output, name => 'output matched expected results');

};

subtest 'SNV One-based test' => sub {
    my $snv_vcf_input = $base_test_dir . "test_snv.vcf";
    my $expected_bed = $base_test_dir . "expected_one_based_snv.bed";
    my $test_output = $temp_base_dir . "test_one_based_snv.bed";

    my $converter = Genome::Model::Tools::Bed::Convert::VcfToBed->create(
       source => $snv_vcf_input,
       output => $test_output,
       one_based => 1,
    );

    ok($converter, 'converter created');
    my $rv = $converter->execute;
    is($rv, 1, 'Testing for succesful execution. Expecting 1. Got : ' . $rv);

    ok(-s $test_output, "output file created");

    compare_ok($expected_bed, $test_output, name => 'output matched expected results');
};

subtest 'No sample data' => sub {
    my $snv_vcf_input = $base_test_dir . "test_no_sample_data.vcf";
    my $expected_bed = $base_test_dir . "expected_no_sample_data.bed";
    my $test_output = $temp_base_dir . "test_no_sample_data.bed";

    my $converter = Genome::Model::Tools::Bed::Convert::VcfToBed->create(
       source => $snv_vcf_input,
       output => $test_output,
       one_based => 1,
    );

    ok($converter, 'converter created');
    my $rv = $converter->execute;
    is($rv, 1, 'Testing for succesful execution. Expecting 1. Got : ' . $rv);

    ok(-s $test_output, "output file created");

    compare_ok($expected_bed, $test_output, name => 'output matched expected results');
};

done_testing();
