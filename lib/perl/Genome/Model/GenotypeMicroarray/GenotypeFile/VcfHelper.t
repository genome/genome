#! /gsc/bin/perl

use strict;
use warnings;

use above 'Genome';

use Test::More;

use_ok('Genome::Model::GenotypeMicroarray::GenotypeFile::VcfHelper') or die;

ok(!eval{Genome::Model::GenotypeMicroarray::GenotypeFile::VcfHelper->header;}, 'Genome::File::Vcf::Header');#, 'header_for_vcf');
isa_ok(Genome::Model::GenotypeMicroarray::GenotypeFile::VcfHelper->header('__TEST_SAMPLE__'), 'Genome::File::Vcf::Header');#, 'header_for_vcf');
is(ref Genome::Model::GenotypeMicroarray::GenotypeFile::VcfHelper->info_order, 'ARRAY', 'info_order');
is(ref Genome::Model::GenotypeMicroarray::GenotypeFile::VcfHelper->info_names_and_ids, 'HASH', 'info_names_and_ids');

done_testing();
