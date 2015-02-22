#!/usr/bin/env genome-perl

use strict;
use warnings;

use above 'Genome';

use Data::Dumper 'Dumper';
use Test::More;

use_ok('Genome::Model::Input') or die;

# CREATE
my $input = Genome::Model::Input->create(
    model_id => 1,
    name => 'fun',
    value_class_name => 'UR::Value',
    value_id => 1,
    filter_desc => 'forward-only',
);
ok($input, 'create input');

# COPY
ok(!$input->copy, 'failed to copy w/o overrides');
like($input->error_message, qr/No overrides to copy model input!/, 'correct error');
ok(!$input->copy(invalid => undef), 'failed to copy w/ invalid overrides');
like($input->error_message, qr/^Unknown overrides given to copy model input! /, 'correct error');
my $copy = $input->copy(model_id => 2);
ok($copy, 'copied input overriding model id'); # part of id
is($copy->model_id, 2, 'correct model_id');
$copy = $input->copy(value_id => 2, filter_desc => undef);
ok($copy, 'copied input overriding value_id and filter id'); # not part of id
is($copy->value_id, 2, 'correct value_id ');
is($copy->filter_desc, undef, 'set filter desc to undef');

done_testing();
