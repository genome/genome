#! /gsc/bin/perl

use strict;
use warnings;

use above 'Genome';

use Test::More;

use_ok('Genome::Model::GenotypeMicroarray::GenotypeFile::Info') or die;

isa_ok(Genome::Model::GenotypeMicroarray::GenotypeFile::Info->header_for_vcf, 'Genome::File::Vcf::Header');#, 'header_for_vcf');
is(ref Genome::Model::GenotypeMicroarray::GenotypeFile::Info->header_for_csv, 'ARRAY', 'header_for_csv');
is(ref Genome::Model::GenotypeMicroarray::GenotypeFile::Info->info_order, 'ARRAY', 'info_order');
is(ref Genome::Model::GenotypeMicroarray::GenotypeFile::Info->info_names_and_ids, 'HASH', 'info_names_and_ids');

done_testing();
