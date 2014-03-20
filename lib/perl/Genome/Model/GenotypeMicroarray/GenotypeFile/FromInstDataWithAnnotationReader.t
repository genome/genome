#! /gsc/bin/perl

use strict;
use warnings;

use above 'Genome';

use Test::More;

use_ok('Genome::Model::GenotypeMicroarray::GenotypeFile::FromInstDataWithAnnotationReader') or die;
use_ok('Genome::File::Vcf::Entry') or die;

sub build_test {
    my ($allele1, $allele2, $reference_allele, $alternate_alleles, $expected_gt, $expected_alts) = @_; 

    my $genotype = { allele1 => $allele1, allele2 => $allele2, };
    my $entry = bless({ reference_allele => $reference_allele, alternate_alleles => $alternate_alleles }, 'Genome::File::Vcf::Entry');

    return sub { 
        my $gt = Genome::Model::GenotypeMicroarray::GenotypeFile::FromInstDataWithAnnotationReader->_gt_for_genotype($genotype, $entry);
        is($gt, $expected_gt, 'correct GT '.$gt.' for alleles: '.$genotype->{allele1}.$genotype->{allele2});
        is_deeply($entry->{alternate_alleles}, $expected_alts, 'correct alternate alleles: '.join('', @{$entry->{alternate_alleles}}));
    };
};

my %params;
subtest 'test dashes' => build_test('-', '-', 'A', [qw/ T /], './.', [qw/ T /]);
subtest 'test homozygous ref' => build_test('A', 'A', 'A', [qw/ T /], '0/0', [qw/ T /]);
subtest 'test homozygous ref w/ many snvs' => build_test('A', 'A', 'A', [qw/ T C G /], '0/0', [qw/ T C G /]);
subtest 'test heterozygous ref' => build_test('A', 'T', 'A', [qw/ T /], '0/1', [qw/ T /]);
subtest 'test heterozygous ref' => build_test('T', 'A', 'A', [qw/ T /], '1/0', [qw/ T /]);
subtest 'test heterozygous ref w/ many snvs' => build_test('T', 'A', 'A', [qw/ C G T /], '3/0', [qw/ C G T /]);
subtest 'test heterozygous ref with new allele' => build_test('A', 'G', 'A', [qw/ T /], '0/2', [qw/ T G /]);

done_testing();
