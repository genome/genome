#!/usr/bin/env genome-perl

use strict;
use warnings;

BEGIN {
    $ENV{UR_DBI_NO_COMMIT} = 1;
};

use Test::More tests => 6;
use above 'Genome';

my $pkg = 'Genome::InstrumentData::Composite::Decorator::Qc';
use_ok($pkg);

my $operation = Genome::WorkflowBuilder::Converge->create(
    name => 'existing operation for test',
    input_properties => ['x'],
    output_properties => [qw(alignment_result y)],
);

my $model = Genome::WorkflowBuilder::DAG->create(
    name => 'model for test',
);

$model->add_operation($operation);

$model->connect_input(
    input_property => 'x',
    destination_property => 'x',
    destination => $operation,
);

$model->connect_output(
    source_property => 'y',
    source => $operation,
    output_property => 'y',
);

my $qc_config = 'test QC config';
my @new_inputs = $pkg->decorate($operation, $model, $qc_config);
is(scalar(@new_inputs), 2, 'got inputs');
is($new_inputs[-1], $qc_config, 'input has correct value');

my $ops = $model->operations;
is(scalar(@$ops), 2, 'an operation was added to the model');
my ($new_op) = grep { $_->isa('Genome::WorkflowBuilder::Command') } @$ops;
is($new_op->command, 'Genome::Qc::Run', 'new operation calls the QC runner');

my @errors = $model->validate;
is(scalar(@errors), 0, 'workflow validates after modifications')
    or diag(join("\n", @errors));
