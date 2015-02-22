#!/usr/bin/env genome-perl

use strict;
use warnings;

use above 'Genome';

use Test::More;

use_ok('Genome::Model::Tools::MetagenomicClassifier::Rdp::ClassifierBase') or die;

class Genome::Model::Tools::MetagenomicClassifier::Rdp::VersionTest {
    is => 'Genome::Model::Tools::MetagenomicClassifier::Rdp::ClassifierBase',
};

# create fails
my $fail = Genome::Model::Tools::MetagenomicClassifier::Rdp::VersionTest->create();
ok(!$fail, 'fail w/o training set');
$fail = Genome::Model::Tools::MetagenomicClassifier::Rdp::VersionTest->create(training_set => 5);
ok(!$fail, 'fail w/ invalid training set');

done_testing();
