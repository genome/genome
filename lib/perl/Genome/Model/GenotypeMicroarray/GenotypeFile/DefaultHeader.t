#! /gsc/bin/perl

use strict;
use warnings;

use above 'Genome';

use Test::More;

use_ok('Genome::Model::GenotypeMicroarray::GenotypeFile::DefaultHeader') or die;

isa_ok(Genome::Model::GenotypeMicroarray::GenotypeFile::DefaultHeader->header, 'Genome::File::Vcf::Header');

done_testing();
