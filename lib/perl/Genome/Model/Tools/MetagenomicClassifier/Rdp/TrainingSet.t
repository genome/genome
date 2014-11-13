#!/usr/bin/env genome-perl

use strict;
use warnings;

use above 'Genome';

require File::Spec;
require Genome::Utility::Test;
use Test::Exception;
use Test::More;

my $class = 'Genome::Model::Tools::MetagenomicClassifier::Rdp::TrainingSet';
use_ok($class) or die;
my $data_dir = Genome::Utility::Test->data_dir_ok('Genome::Model::Tools::MetagenomicClassifier::Rdp');
my $training_set_path = File::Spec->catfile($data_dir, 'training-set');

my $set_name = 10;
my @valid_set_names = $class->valid_set_names;
ok(@valid_set_names, 'valid set names');
my $path_for_set_name = $class->path_for_set_name($valid_set_names[0]);
is($path_for_set_name, File::Spec->catfile($class->base_path, $valid_set_names[0]), 'path for set name');

my $ts = $class->create(path => $training_set_path);
ok($ts, "create training set");
ok($ts->classifier_properties_path, 'classifier properties path');

# Failures
$ts = $class->create(path => '0');
ok($ts, "create w/ invalid training path");
throws_ok(
    sub{ $ts->classifier_properties_path; },
    qr/No classifier properties \(rRNAClassifier.properties\) in training path! 0/,
    'Correct error for non existing classifier_properties_path',
);
throws_ok(
    sub{ $ts->taxonomy_path; },
    qr/No taxonomy XML \(bergeyTrainingTree.xml\) in training path! 0/,
    'Correct error for non existing taxonomy_path',
);

throws_ok(
    sub{ $class->path_for_set_name() },
    qr/No set name given to get path!/,
    'Failed to get path w/o set name',
);

throws_ok(
    sub{ $class->path_for_set_name(0) },
    qr/Invalid set name! 0/,
    'Failed to get path w/ invalid set name',
);

done_testing();
