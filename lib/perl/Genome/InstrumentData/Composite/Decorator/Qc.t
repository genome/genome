#!/usr/bin/env genome-perl

use strict;
use warnings;

BEGIN {
    $ENV{UR_DBI_NO_COMMIT} = 1;
};

use Test::More tests => 6;
use above 'Genome';
use Workflow;

my $pkg = 'Genome::InstrumentData::Composite::Decorator::Qc';
use_ok($pkg);

my $operation = Workflow::Operation->create(
    name => 'existing operation for test',
    operation_type => Workflow::OperationType::Converge->create(
        input_properties => ['x'],
        output_properties => [qw(alignment_result y)],
    ),
);

my $model = Workflow::Model->create(
    name => 'model for test',
    input_properties => [qw(result_users x)],
    optional_input_properties => [],
    output_properties => ['y']
);

$operation->workflow_model($model);

$model->add_link(
    left_property => 'x',
    left_operation => $model->get_input_connector,
    right_property => 'x',
    right_operation => $operation,
);

$model->add_link(
    left_property => 'y',
    left_operation => $operation,
    right_property => 'y',
    right_operation => $model->get_output_connector,
);

my $qc_config = 'test QC config';
my @new_inputs = $pkg->decorate($operation, $qc_config);
is(scalar(@new_inputs), 2, 'got inputs');
is($new_inputs[-1], $qc_config, 'input has correct value');

my @ops = $model->operations;
is(scalar(@ops), 4, 'an operation was added to the model');
my ($new_op) = grep { $_->operation_type->isa('Workflow::OperationType::Command') } @ops;
is($new_op->operation_type->command_class_name, 'Genome::Qc::Run', 'new operation calls the QC runner');

my @errors = $model->validate;
is(scalar(@errors), 0, 'workflow validates after modifications')
    or diag(join("\n", @errors));
