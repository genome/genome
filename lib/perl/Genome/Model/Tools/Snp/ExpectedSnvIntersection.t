#!/usr/bin/env genome-perl

use strict;
use warnings;

use Test::More skip_all => 'in development'; #tests => 5;
use File::Slurp;

use above 'Genome';

BEGIN {
        use_ok('Genome::Model::Tools::Snp::ExpectedSnvIntersection');
    };

my $dir = $ENV{GENOME_TEST_INPUTS} . '/Genome-Model-Tools-Snp/ExpectedSnvIntersection';

my $expected_snv_file = "$dir/expected.snv";
my $maq_snv_file  = "$dir/maq.snv";
my $sam_snv_file  = "$dir/sam.snv";

my $exp_maq_report = read_file($dir.'/report.maq.ori');
my $exp_sam_report = read_file($dir.'/report.sam.ori');

=cut

my $maq_gsi = Genome::Model::Tools::Snp::ExpectedSnvIntersection->create(
    snv_file      => $maq_snv_file,
    expected_snv_file => $expected_snv_file,
);

isa_ok($maq_gsi,'Genome::Model::Tools::Snp::ExpectedSnvIntersection');

=cut

my $maq_report = `gt snv expected-snv-intersection --snv-file $maq_snv_file --expected-snv-file $expected_snv_file`;

ok($maq_report,'MAQ expected-snv-intersection execute ok');
is($maq_report, $exp_maq_report, 'MAQ expected-snv-intersection output matches the expected original one.');

=cut

my $sam_gsi = Genome::Model::Tools::Snp::ExpectedSnvIntersection->create(
    snv_file      => $sam_snv_file,
    expected_snv_file => $expected_snv_file,
    sam_format    => 1,
);

my $sam_report = $sam_gsi->execute;

=cut

my $sam_report = `gt snv expected-snv-intersection --snv-file $sam_snv_file --expected-snv-file $expected_snv_file --snv-format sam`;

ok($sam_report,'SAM expected-snv-intersection execute ok');
is($sam_report, $exp_sam_report, 'SAM expected-snv-intersection output matches the expected original one.');
