#!/usr/bin/env genome-perl

use strict;
use warnings;

use above 'Genome';

use File::Compare;
use Test::More tests => 7;

BEGIN {
    use_ok('Genome::Model::Tools::Sam::LimitVariants');
}

my $data_dir = $ENV{GENOME_TEST_INPUTS} . '/Genome-Model-Tools-Sam-LimitVariants';
my $tmp_dir  = Genome::Sys->create_temp_directory('Genome-Model-Tools-Sam-LimitVariants-'. Genome::Sys->username);

my $bed_file = $data_dir .'/test_regions.bed';
my $snp_file = $data_dir .'/test.snp';
my $expected_snp_file = $data_dir .'/expected.snp';
my $output_snp_file = $tmp_dir .'/output.snp';

my $limit_snps = Genome::Model::Tools::Sam::LimitVariants->create(
    variants_file   => $snp_file,
    output_file   => $output_snp_file,
    bed_file => $bed_file,
);
isa_ok($limit_snps,'Genome::Model::Tools::Sam::LimitVariants');
ok($limit_snps->execute,'executed ok');
cmp_ok(compare($output_snp_file, $expected_snp_file), '==', 0, 'LimitVaraints file was created ok');

my $indel_file = $data_dir .'/test.indel';
my $expected_indel_file = $data_dir .'/expected.indel';
my $output_indel_file = $tmp_dir .'/output.indel';
my $limit_indels = Genome::Model::Tools::Sam::LimitVariants->create(
    variants_file   => $indel_file,
    output_file   => $output_indel_file,
    bed_file => $bed_file,
);
isa_ok($limit_indels,'Genome::Model::Tools::Sam::LimitVariants');
ok($limit_indels->execute,'executed ok');
cmp_ok(compare($output_indel_file, $expected_indel_file), '==', 0, 'LimitVaraints file was created ok');




exit;

