#!/usr/bin/env genome-perl

use strict;
use warnings;

use Test::More;

if (Genome::Config->arch_os ne 'x86_64') {
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

subtest 'Indel simplification - none' => sub {
    my $ref = 'G';
    my $var = undef;

    my ($result_ref, $result_var, $shift) = $instance->_simplify_indel_allele($ref, $var);
    is($result_ref, 'G', "ref allele as expected");
    is($result_var, '', "var allele as expected");
    is($shift, 0, "no position shift as expected");
};

subtest 'Indel simplification - simple insertion' => sub {
    my $ref = 'G';
    my $var = 'GT';

    my ($result_ref, $result_var, $shift) = $instance->_simplify_indel_allele($ref, $var);
    is($result_ref, q{}, "ref allele as expected");
    is($result_var, 'T', "var allele as expected");
    is($shift, 1, "position shift as expected");
};

subtest 'Indel simplification - simple deletion' => sub {
    my $ref = 'GT';
    my $var = 'G';

    my ($result_ref, $result_var, $shift) = $instance->_simplify_indel_allele($ref, $var);
    is($result_ref, 'T', "ref allele as expected");
    is($result_var, "", "var allele as expected");
    is($shift, 1, "position shift as expected");
};

subtest 'Indel simplification - ambiguous deletion' => sub {
    my $ref = 'GTTC';
    my $var = 'GTC';

    my ($result_ref, $result_var, $shift) = $instance->_simplify_indel_allele($ref, $var);
    is($result_ref, 'T', "ref allele as expected");
    is($result_var, "", "var allele as expected");
    is($shift, 1, "position shift as expected");
};

subtest 'Indel simplification - ambiguous insertion' => sub {
    my $ref = 'AGC';
    my $var = 'AGGGGGC';

    my ($result_ref, $result_var, $shift) = $instance->_simplify_indel_allele($ref, $var);
    is($result_ref, "", "ref allele as expected");
    is($result_var, "GGGG", "var allele as expected");
    is($shift, 1, "position shift as expected");
};

subtest 'Indel simplification - ambiguous deletion' => sub {
    my $ref = 'AGTGTC';
    my $var = 'AGTC';

    my ($result_ref, $result_var, $shift) = $instance->_simplify_indel_allele($ref, $var);
    is($result_ref, "GT", "ref allele as expected");
    is($result_var, "", "var allele as expected");
    is($shift, 1, "position shift as expected");
};

subtest 'Indel gt - simple insertion' => sub {
    my $ref = 'G';
    my @var = ('GT');

    my ($result, $shifts) = $instance->_convert_indel_gt_to_bed($ref, @var);
    is_deeply($result, [['*', 'T']], "alleles as expected");
    is_deeply($shifts, [1], "no position shift as expected");
};

subtest 'Indel gt - simple deletion' => sub {
    my $ref = 'GT';
    my @var = ('G');

    my ($result, $shifts) = $instance->_convert_indel_gt_to_bed($ref, @var);
    is_deeply($result, [['T', '*']], "alleles as expected");
    is_deeply($shifts, [1], "no position shift as expected");
};

subtest 'Indel gt - ambiguous deletion' => sub {
    my $ref = 'GTTTC';
    my @var = ('GTC');

    my ($result, $shifts) = $instance->_convert_indel_gt_to_bed($ref, @var);
    is_deeply($result, [['TT', '*']], "alleles as expected");
    is_deeply($shifts, [1], "no position shift as expected");
};

subtest 'Indel gt - ambiguous insertion' => sub {
    my $ref = 'AGC';
    my @var = ('AGGGGGC');

    my ($result, $shifts) = $instance->_convert_indel_gt_to_bed($ref, @var);
    is_deeply($result, [['*', 'GGGG']], "alleles as expected");
    is_deeply($shifts, [1], "no position shift as expected");
};

subtest 'Indel gt - ambiguous insertion' => sub {
    my $ref = 'AGTC';
    my @var = ('AGTGTC');

    # GT--C
    # GTGTC

    my ($result, $shifts) = $instance->_convert_indel_gt_to_bed($ref, @var);
    is_deeply($result, [['*', 'GT']], "alleles as expected");
    is_deeply($shifts, [1], "position shift as expected");
};

subtest 'Indel gt - ambiguous deletion' => sub {
    my $ref = 'AGTGTC';
    my @var = ('AGTC');

    # GT--C
    # GTGTC

    my ($result, $shifts) = $instance->_convert_indel_gt_to_bed($ref, @var);
    is_deeply($result, [['GT', '*']], "alleles as expected");
    is_deeply($shifts, [1], "no position shift as expected");
};

my $base_test_dir = $ENV{GENOME_TEST_INPUTS} . "/Genome-Model-Tools-Bed-Convert-Vcf-To-Bed/v1/";
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

done_testing();
