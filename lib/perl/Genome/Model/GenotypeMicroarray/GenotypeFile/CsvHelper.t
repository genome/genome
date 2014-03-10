#! /gsc/bin/perl

use strict;
use warnings;

use above 'Genome';

use Test::More;

use_ok('Genome::Model::GenotypeMicroarray::GenotypeFile::CsvHelper') or die;

is(ref Genome::Model::GenotypeMicroarray::GenotypeFile::CsvHelper->column_names, 'ARRAY', 'colum names');

done_testing();
