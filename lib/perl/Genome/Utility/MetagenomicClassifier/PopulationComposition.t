#!/usr/bin/env genome-perl

use strict;
use warnings;

use above 'Genome';

use Storable qw/ retrieve /;
use Test::More;

use_ok('Genome::Utility::MetagenomicClassifier::PopulationComposition') or die;
use_ok('Genome::Utility::MetagenomicClassifier::SequenceClassification') or die;

my $population_composition = Genome::Utility::MetagenomicClassifier::PopulationComposition->new( 
    confidence_threshold => .8,
);
ok($population_composition, 'created population composition');

can_ok($population_composition, 'add_classification');
my $classifications = retrieve($ENV{GENOME_TEST_INPUTS} . '/Genome-Utility-MetagenomicClassifier/classifications.stor');
for my $classification ( @$classifications ) { # should be 10
    $population_composition->add_classification($classification)
}

my $eval = eval {
    Genome::Utility::MetagenomicClassifier::PopulationComposition->new(confidence_threshold => '1.5.5');
};
diag("$@\n");
ok(!$eval, 'Failed as expected - create w/ invalid confidence_threhold');

done_testing();
