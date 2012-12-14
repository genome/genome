#!/usr/bin/env genome-perl

use strict;
use warnings;

use above 'Genome';

use Test::More;

use_ok('Genome::Utility::MetagenomicClassifier::Rdp::Version2x1') or die;
use_ok('Genome::Utility::MetagenomicClassifier::PopulationCompositionFactory') or die;
use_ok('Genome::Utility::MetagenomicClassifier::Rdp') or die;

my $classifier = Genome::Utility::MetagenomicClassifier::Rdp::Version2x1->new(
    training_set => 'broad',
);
ok($classifier, 'Created rdp classifier');
my $factory = Genome::Utility::MetagenomicClassifier::PopulationCompositionFactory->instance;
ok($factory, 'Got factory instance');
my $composition = $factory->get_composition(
    classifier => $classifier,
    fasta_file => $ENV{GENOME_TEST_INPUTS} . '/Genome-Utility-MetagenomicClassifier/U_PR-JP_TS1_2PCA.fasta',
);
ok($composition, 'Got composition from factory');
isa_ok($composition, 'Genome::Utility::MetagenomicClassifier::PopulationComposition');

done_testing();
