#!/usr/bin/env genome-perl

use strict;
use warnings;

use above 'Genome';

require File::Spec;
require Genome::Utility::Test;
use Test::More;

my $class = 'Genome::Model::Tools::MetagenomicClassifier::Rdp::ListGenera';
use_ok($class) or die;
my $data_dir = Genome::Utility::Test->data_dir_ok('Genome::Model::Tools::MetagenomicClassifier::Rdp');
my $training_set_path = File::Spec->catfile($data_dir, 'training-set');

my $list = $class->create(training_set_path => $training_set_path);
ok($list, 'create w/ traning set path');
ok($list->execute, 'execute w/ traning set path');

$list = $class->create(training_set_name => '10');
ok($list, 'create w/ training set name');

# FAILS
$list = $class->create();
ok(!$list->execute, 'failed to execute w/o params');

$list = $class->create(training_set_path => 0);
ok(!$list->execute, 'failed to execute w/ invalid training set path');

done_testing();
