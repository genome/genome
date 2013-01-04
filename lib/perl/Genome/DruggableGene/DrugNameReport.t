#!/usr/bin/env genome-perl

use strict;
use warnings;

use above 'Genome';
use Test::More tests => 15;

SKIP: {
skip "DB in flux", 15;
use_ok('Genome::FeatureList');

my $drug = Genome::DruggableGene::DrugNameReport->get('0007e48149af4a19a08fba1dfa02f8f9');
ok($drug->citation, 'Drug has a citation');
ok($drug->original_data_source_url, 'Drug has an orignal data source url');

#withdrawn
my $withdrawn = Genome::DruggableGene::DrugNameReport->get('01269b0484cb4bd8bc20a403c70840b9');
ok($withdrawn, 'got a withdrawn drug');
ok($withdrawn->is_withdrawn, 'drug is withdrawn');
my $not_withdrawn = Genome::DruggableGene::DrugNameReport->get('0007e481-49af-4a19-a08f-ba1dfa02f8f9');
ok($not_withdrawn, 'got a non withdrawn drug');
ok(!$not_withdrawn->is_withdrawn, 'drug is not withdrawn');

#nutraceutical
my $nutra = Genome::DruggableGene::DrugNameReport->get('0252b667d40d4e00852a782bea58bd9a');
ok($nutra, 'got a drug nutraceutical');
ok($nutra->is_nutraceutical, 'drug is nutraceutical');

my $not_nutraceutical = Genome::DruggableGene::DrugNameReport->get('0007e481-49af-4a19-a08f-ba1dfa02f8f9');
ok($not_nutraceutical, 'got a non nutraceutical drug');
ok(!$not_nutraceutical->is_nutraceutical, 'drug is not nutraceutical');

#approved
my $approved = Genome::DruggableGene::DrugNameReport->get('0047db3602ce4fcb80d18ca87fa671bc');
ok($approved, 'got a drug approved');
ok($approved->is_approved, 'drug is approved');

my $not_approved = Genome::DruggableGene::DrugNameReport->get('0007e481-49af-4a19-a08f-ba1dfa02f8f9');
ok($not_approved, 'got a non approved drug');
ok(!$not_approved->is_approved, 'drug is not nutraceutical');
}

1;
