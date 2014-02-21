#! /gsc/bin/perl

use strict;
use warnings;

use above 'Genome';

use Test::More;

use_ok('Genome::Model::GenotypeMicroarray::GenotypeFile::DefaultHeader') or die;

is(ref Genome::Model::GenotypeMicroarray::GenotypeFile::DefaultHeader->header_lines, 'ARRAY', 'header_lines');
isa_ok(Genome::Model::GenotypeMicroarray::GenotypeFile::DefaultHeader->header, 'Genome::File::Vcf::Header');
is(ref Genome::Model::GenotypeMicroarray::GenotypeFile::DefaultHeader->info_order, 'ARRAY', 'info_order');
is(ref Genome::Model::GenotypeMicroarray::GenotypeFile::DefaultHeader->info_names_and_ids, 'HASH', 'info_names_and_ids');

done_testing();
