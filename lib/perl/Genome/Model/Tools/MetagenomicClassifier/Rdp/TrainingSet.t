#!/usr/bin/env genome-perl

use strict;
use warnings;

use above 'Genome';
use Test::Exception;
use Test::More;

my $class = 'Genome::Model::Tools::MetagenomicClassifier::Rdp::TrainingSet';
use_ok($class) or die;

my $training_set_name = 10;
my $ts_10 = $class->create_from_training_set_name($training_set_name);
ok($ts_10, "create training set $training_set_name");
my $training_path = $ts_10->training_path;
ok($training_path, "path for training set $training_set_name");
ok($ts_10->classifier_properties_path, 'classifier properties path');

# Failures
throws_ok(
    sub{ $class->create_from_training_set_name() },
    qr/No training set given to get training path!/,
    'Failed to create without training set name',
);

throws_ok(
    sub{ $class->create_from_training_set_name(0) },
    qr/Invalid training set \(0\) given to create_from_training_set_name!/,
    'Failed to create without training set name',
);

done_testing();
