use strict;
use warnings;

use above 'Genome';
use Test::More;
use Test::Exception;


use_ok('Genome::WorkflowBuilder::DAG');

subtest 'Simple DAG' => sub {
    my $dag = Genome::WorkflowBuilder::DAG->create(
        name => 'top level',
        log_dir => '/tmp',
    );

    my $op = Genome::WorkflowBuilder::Command->create(
        name => 'some op',
        command => 'Genome::WorkflowBuilder::Test::DummyCommand',
    );
    $dag->add_operation($op);

    $dag->connect_input(
        input_property => 'some_external_input',
        destination => $op,
        destination_property => 'input',
    );
    $dag->connect_output(
        output_property => 'some_external_output',
        source => $op,
        source_property => 'single_output',
    );

    my $expected_xml = <<EOS;
<?xml version="1.0"?>
<operation name="top level" logDir="/tmp">
  <operationtype typeClass="Workflow::OperationType::Model">
    <inputproperty>some_external_input</inputproperty>
    <outputproperty>some_external_output</outputproperty>
  </operationtype>
  <operation name="some op">
    <operationtype typeClass="Workflow::OperationType::Command" lsfQueue="$ENV{GENOME_LSF_QUEUE_BUILD_WORKER_ALT}" lsfResource="-M 25000000 -R 'select[mem&gt;25000] rusage[mem=25000]'" commandClass="Genome::WorkflowBuilder::Test::DummyCommand">
      <inputproperty>input</inputproperty>
      <outputproperty>many_output</outputproperty>
      <outputproperty>result</outputproperty>
      <outputproperty>single_output</outputproperty>
    </operationtype>
  </operation>
  <link fromOperation="input connector" fromProperty="some_external_input" toOperation="some op" toProperty="input"/>
  <link fromOperation="some op" fromProperty="single_output" toOperation="output connector" toProperty="some_external_output"/>
</operation>
EOS

    cmp_xml($dag->get_xml, $expected_xml);
};

subtest 'Invalid DAG Name' => sub {
    my $dag = Genome::WorkflowBuilder::DAG->create(name => 'input connector');

    eval {
        diag "Expect one error message about invalid operation name:";
        $dag->validate;
    };

    ok($@, 'invalid operation name fails to validate');
};

subtest 'Non-Unique Operation Names' => sub {
    my $dag = Genome::WorkflowBuilder::DAG->create(name => 'top level');

    my $op1 = Genome::WorkflowBuilder::Command->create(
        name => 'duplicate name',
        command => 'Genome::WorkflowBuilder::Test::DummyCommand'
    );
    $dag->add_operation($op1);

    my $op2 = Genome::WorkflowBuilder::Command->create(
        name => 'duplicate name',
        command => 'Genome::WorkflowBuilder::Test::DummyCommand'
    );
    $dag->add_operation($op2);

    eval {
        diag "Expect one error message about duplicate operation name";
        $dag->validate;
    };

    ok($@, 'duplicate operation names not allowed');
};

subtest 'Unowned Operations' => sub {
    my $dag = Genome::WorkflowBuilder::DAG->create(name => 'top level');

    my $owned_op = Genome::WorkflowBuilder::Command->create(
        name => 'owned op',
        command => 'Genome::WorkflowBuilder::Test::DummyCommand'
    );
    $dag->add_operation($owned_op);

    my $unowned_op = Genome::WorkflowBuilder::Command->create(
        name => 'unowned op',
        command => 'Genome::WorkflowBuilder::Test::DummyCommand'
    );

    $dag->create_link(
        source => $owned_op,
        source_property => 'single_output',
        destination => $unowned_op,
        destination_property => 'input'
    );

    eval {
        diag "Expect one error message about unowned operations";
        $dag->validate;
    };

    ok($@, 'unowned operations not allowed in dag links');
};

subtest 'Mandatory Inputs' => sub {
    my $dag = Genome::WorkflowBuilder::DAG->create(name => 'top level');

    my $op = Genome::WorkflowBuilder::Command->create(
        name => 'some op',
        command => 'Genome::WorkflowBuilder::Test::DummyCommand'
    );
    $dag->add_operation($op);

    eval {
        diag 'Expect one error message about missing operation inputs';
        $dag->validate();
    };

    ok($@, 'missing mandatory inputs not allowed in dag');
};

subtest 'Conflicting Inputs' => sub {
    my $dag = Genome::WorkflowBuilder::DAG->create(name => 'top level');

    my $op = Genome::WorkflowBuilder::Command->create(
        name => 'some op',
        command => 'Genome::WorkflowBuilder::Test::DummyCommand'
    );
    $dag->add_operation($op);

    $dag->connect_input(
        input_property => 'external_input_a',
        destination => $op,
        destination_property => 'input',
    );

    $dag->connect_input(
        input_property => 'external_input_b',
        destination => $op,
        destination_property => 'input',
    );

    eval {
        diag 'Expect one error message about conflicting operation inputs';
        $dag->validate();
    };

    ok($@, 'conflicting inputs not allowed in dag');
};

subtest 'XML Round Trip' => sub {
    my $xml = <<EOS;
<?xml version="1.0"?>
<operation name="top level" parallelBy="some_external_input" logDir="/tmp">
  <operationtype typeClass="Workflow::OperationType::Model">
    <inputproperty>some_external_input</inputproperty>
    <outputproperty>some_external_output</outputproperty>
  </operationtype>
  <operation name="some op">
    <operationtype typeClass="Workflow::OperationType::Command" lsfQueue="$ENV{GENOME_LSF_QUEUE_BUILD_WORKER_ALT}" lsfResource="-M 25000000 -R 'select[mem&gt;25000] rusage[mem=25000]'" commandClass="Genome::WorkflowBuilder::Test::DummyCommand">
      <inputproperty>input</inputproperty>
      <outputproperty>many_output</outputproperty>
      <outputproperty>result</outputproperty>
      <outputproperty>single_output</outputproperty>
    </operationtype>
  </operation>
  <link fromOperation="input connector" fromProperty="some_external_input" toOperation="some op" toProperty="input"/>
  <link fromOperation="some op" fromProperty="single_output" toOperation="output connector" toProperty="some_external_output"/>
</operation>
EOS

    my $dag = Genome::WorkflowBuilder::DAG->from_xml($xml);
    cmp_xml($dag->get_xml, $xml);
};

subtest 'Converge XML Round Trip' => sub {
    my $xml = <<EOS;
<?xml version="1.0"?>
<operation name="top level" parallelBy="some_external_input" logDir="/tmp">
  <operationtype typeClass="Workflow::OperationType::Model">
    <inputproperty>external_input_0</inputproperty>
    <inputproperty>external_input_1</inputproperty>
    <outputproperty>external_output</outputproperty>
  </operationtype>
  <operation name="some op">
    <operationtype typeClass="Workflow::OperationType::Converge">
      <inputproperty>input_0</inputproperty>
      <inputproperty>input_1</inputproperty>
      <outputproperty>converge_output</outputproperty>
      <outputproperty>result</outputproperty>
    </operationtype>
  </operation>
  <link fromOperation="input connector" fromProperty="external_input_0" toOperation="some op" toProperty="input_0"/>
  <link fromOperation="input connector" fromProperty="external_input_1" toOperation="some op" toProperty="input_1"/>
  <link fromOperation="some op" fromProperty="converge_output" toOperation="output connector" toProperty="external_output"/>
</operation>
EOS

    my $dag = Genome::WorkflowBuilder::DAG->from_xml($xml);
    cmp_xml($dag->get_xml, $xml);
};

subtest 'Converge DAG' => sub {
    my $xml = <<EOS;
<?xml version="1.0"?>
<operation name="top level">
  <operationtype typeClass="Workflow::OperationType::Model">
    <inputproperty>external_input_0</inputproperty>
    <inputproperty>external_input_1</inputproperty>
    <outputproperty>external_output</outputproperty>
  </operationtype>
  <operation name="some op">
    <operationtype typeClass="Workflow::OperationType::Converge">
      <inputproperty>input_1</inputproperty>
      <inputproperty>input_0</inputproperty>
      <outputproperty>converge_output</outputproperty>
    </operationtype>
  </operation>
  <link fromOperation="input connector" fromProperty="external_input_0" toOperation="some op" toProperty="input_0"/>
  <link fromOperation="input connector" fromProperty="external_input_1" toOperation="some op" toProperty="input_1"/>
  <link fromOperation="some op" fromProperty="converge_output" toOperation="output connector" toProperty="external_output"/>
</operation>
EOS

    my $dag = Genome::WorkflowBuilder::DAG->create(name => 'top level');
    my $converge = $dag->add_operation(
        Genome::WorkflowBuilder::Converge->create(name => 'some op'));
    $dag->connect_input(input_property => 'external_input_1',
        destination => $converge, destination_property => 'input_1');
    $dag->connect_input(input_property => 'external_input_0',
        destination => $converge, destination_property => 'input_0');

    $dag->connect_output(output_property => 'external_output',
        source => $converge, source_property => 'converge_output');

    cmp_xml($dag->get_xml, $xml);
};

subtest 'Nested DAG with constant input' => sub {
    my $outer = Genome::WorkflowBuilder::DAG->create(
        name => 'outer',
        log_dir => '/tmp',
    );

    my $inner = Genome::WorkflowBuilder::DAG->create(
        name => 'inner',
        log_dir => '/tmp',
    );

    my $op = Genome::WorkflowBuilder::Command->create(
        name => 'some op',
        command => 'Genome::WorkflowBuilder::Test::DummyCommand',
    );

    dies_ok {$op->declare_constant('non-property' => 'value')}
        'cannot declare constants that are not input properties';
    $op->declare_constant(input => 'foo');
    lives_ok {$op->declare_constant(input => 'bar')}
        'can declare constants more than once';

    $inner->add_operation($op);

    # no 'connect_input' needed
    $inner->connect_output(
        output_property => 'output',
        source => $op,
        source_property => 'single_output',
    );

    $outer->add_operation($inner);
    # no 'connect_input' needed
    $outer->connect_output(
        output_property => 'output',
        source => $inner,
        source_property => 'output',
    );

    my $expected_xml = <<EOS;
<?xml version="1.0"?>
<operation name="outer" logDir="/tmp">
  <operationtype typeClass="Workflow::OperationType::Model">
    <inputproperty>inner.some op.input</inputproperty>
    <outputproperty>output</outputproperty>
  </operationtype>
  <operation name="inner" logDir="/tmp">
    <operationtype typeClass="Workflow::OperationType::Model">
      <inputproperty>some op.input</inputproperty>
      <outputproperty>output</outputproperty>
    </operationtype>
    <operation name="some op">
      <operationtype typeClass="Workflow::OperationType::Command" lsfQueue="apipe" lsfResource="-M 25000000 -R 'select[mem&gt;25000] rusage[mem=25000]'" commandClass="Genome::WorkflowBuilder::Test::DummyCommand">
        <inputproperty>input</inputproperty>
        <outputproperty>many_output</outputproperty>
        <outputproperty>result</outputproperty>
        <outputproperty>single_output</outputproperty>
      </operationtype>
    </operation>
    <link fromOperation="input connector" fromProperty="some op.input" toOperation="some op" toProperty="input"/>
    <link fromOperation="some op" fromProperty="single_output" toOperation="output connector" toProperty="output"/>
  </operation>
  <link fromOperation="inner" fromProperty="output" toOperation="output connector" toProperty="output"/>
  <link fromOperation="input connector" fromProperty="inner.some op.input" toOperation="inner" toProperty="some op.input"/>
</operation>
EOS
    is_deeply($outer->constant_values, {'inner.some op.input' => 'bar'},
        'found expected constants');
    cmp_xml($outer->get_xml, $expected_xml);
};


done_testing();

sub cmp_xml {
    my ($got, $expected) = @_;

    is_deeply([split(/\n/, $got)], [split(/\n/, $expected)], 'got expected xml');
}
