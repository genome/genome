#!/usr/bin/env genome-perl

use strict;
use warnings 'FATAL';

use above 'Genome';

use Genome::Test::Factory::Sample;
use Test::More;
plan tests => 4;

my $class = 'Genome::Model::ReferenceAlignment::Command::GenotypeMicroarrayConcordance';
use_ok($class) or die;

my $sample = Genome::Test::Factory::Sample->generate_obj(name => '_TEST(1)_SAMPLE_');

my $gc = $class->create();
ok($gc, 'crearte genotype concordance command');

my $sorted_vcf = $gc->sorted_microarray_vcf_for_genotype_sample($sample);
my $sanitized_sample_name = Genome::Utility::Text::sanitize_string_for_filesystem($sample->name);
like($sorted_vcf, qr/$sanitized_sample_name/, 'sorted_microarray_vcf_for_genotype_sample has sanitized sample name');
is($sorted_vcf, $gc->sorted_microarray_vcf_for_genotype_sample($sample), 'same sorted_microarray_vcf_for_genotype_sample retrieved again');

done_testing();
