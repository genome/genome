#!/usr/bin/env genome-perl

use strict;
use warnings;
use File::Slurp "read_file";
use Test::More;

use above 'Genome';
use Genome::Utility::Vcf qw(get_vcf_header convert_file_with_alpha_gt_values_to_numeric diff_vcf_file_vs_file convert_indel_gt_to_bed);
                    

# Test header reading
{
    my $input_filename = join( '/',  __FILE__ . ".d",  'input.clean.vcf');
    my $expected_filename = join( '/',  __FILE__ . ".d",  'expected.txt');

    my $header = get_vcf_header($input_filename);
    my $expected_header = read_file($expected_filename);

    is($header, $expected_header, "Read in the header of a vcf file correctly.");
}

# Test alpha GT value -> numeric alt value fixing from polymutt denovo files
{
    my $input_filename = join( '/',  __FILE__ . ".d",  'input_denovo.vcf.gz');
    my $expected_filename = join( '/',  __FILE__ . ".d",  'expected_denovo.vcf.gz');

    my $output_filename = Genome::Sys->create_temp_file_path;
    convert_file_with_alpha_gt_values_to_numeric($input_filename, $output_filename);

    my $diff = diff_vcf_file_vs_file($expected_filename, $output_filename);
    is($diff, "", "convert_file_with_alpha_gt_values_to_numeric outputs as expected");
}

subtest 'Indel simplification - none' => sub {
    my $ref = 'G';
    my $var = undef;

    my ($result_ref, $result_var, $shift) = Genome::Utility::Vcf::_simplify_indel_allele($ref, $var);
    is($result_ref, 'G', "ref allele as expected");
    is($result_var, '', "var allele as expected");
    is($shift, 0, "no position shift as expected");
};

subtest 'Indel simplification - simple insertion' => sub {
    my $ref = 'G';
    my $var = 'GT';

    my ($result_ref, $result_var, $shift) = Genome::Utility::Vcf::_simplify_indel_allele($ref, $var);
    is($result_ref, q{}, "ref allele as expected");
    is($result_var, 'T', "var allele as expected");
    is($shift, 1, "position shift as expected");
};

subtest 'Indel simplification - simple deletion' => sub {
    my $ref = 'GT';
    my $var = 'G';

    my ($result_ref, $result_var, $shift) = Genome::Utility::Vcf::_simplify_indel_allele($ref, $var);
    is($result_ref, 'T', "ref allele as expected");
    is($result_var, "", "var allele as expected");
    is($shift, 1, "position shift as expected");
};

subtest 'Indel simplification - ambiguous deletion' => sub {
    my $ref = 'GTTC';
    my $var = 'GTC';

    my ($result_ref, $result_var, $shift) = Genome::Utility::Vcf::_simplify_indel_allele($ref, $var);
    is($result_ref, 'T', "ref allele as expected");
    is($result_var, "", "var allele as expected");
    is($shift, 1, "position shift as expected");
};

subtest 'Indel simplification - ambiguous insertion' => sub {
    my $ref = 'AGC';
    my $var = 'AGGGGGC';

    my ($result_ref, $result_var, $shift) = Genome::Utility::Vcf::_simplify_indel_allele($ref, $var);
    is($result_ref, "", "ref allele as expected");
    is($result_var, "GGGG", "var allele as expected");
    is($shift, 1, "position shift as expected");
};

subtest 'Indel simplification - ambiguous deletion' => sub {
    my $ref = 'AGTGTC';
    my $var = 'AGTC';

    my ($result_ref, $result_var, $shift) = Genome::Utility::Vcf::_simplify_indel_allele($ref, $var);
    is($result_ref, "GT", "ref allele as expected");
    is($result_var, "", "var allele as expected");
    is($shift, 1, "position shift as expected");
};

subtest 'Indel gt - simple insertion' => sub {
    my $ref = 'G';
    my @var = ('GT');

    my ($result, $shifts) = convert_indel_gt_to_bed($ref, @var);
    is_deeply($result, [['*', 'T']], "alleles as expected");
    is_deeply($shifts, [1], "no position shift as expected");
};

subtest 'Indel gt - simple deletion' => sub {
    my $ref = 'GT';
    my @var = ('G');

    my ($result, $shifts) = convert_indel_gt_to_bed($ref, @var);
    is_deeply($result, [['T', '*']], "alleles as expected");
    is_deeply($shifts, [1], "no position shift as expected");
};

subtest 'Indel gt - ambiguous deletion' => sub {
    my $ref = 'GTTTC';
    my @var = ('GTC');

    my ($result, $shifts) = convert_indel_gt_to_bed($ref, @var);
    is_deeply($result, [['TT', '*']], "alleles as expected");
    is_deeply($shifts, [1], "no position shift as expected");
};

subtest 'Indel gt - ambiguous insertion' => sub {
    my $ref = 'AGC';
    my @var = ('AGGGGGC');

    my ($result, $shifts) = convert_indel_gt_to_bed($ref, @var);
    is_deeply($result, [['*', 'GGGG']], "alleles as expected");
    is_deeply($shifts, [1], "no position shift as expected");
};

subtest 'Indel gt - ambiguous insertion' => sub {
    my $ref = 'AGTC';
    my @var = ('AGTGTC');

    # GT--C
    # GTGTC

    my ($result, $shifts) = convert_indel_gt_to_bed($ref, @var);
    is_deeply($result, [['*', 'GT']], "alleles as expected");
    is_deeply($shifts, [1], "position shift as expected");
};

subtest 'Indel gt - ambiguous deletion' => sub {
    my $ref = 'AGTGTC';
    my @var = ('AGTC');

    # GT--C
    # GTGTC

    my ($result, $shifts) = convert_indel_gt_to_bed($ref, @var);
    is_deeply($result, [['GT', '*']], "alleles as expected");
    is_deeply($shifts, [1], "no position shift as expected");
};
done_testing();
1;
