use strict;
use warnings;

use above 'Genome';
use Test::More;


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
    <operationtype typeClass="Workflow::OperationType::Command" lsfQueue="apipe" lsfResource="-M 25000000 -R 'select[mem&gt;25000] rusage[mem=25000]'" commandClass="Genome::WorkflowBuilder::Test::DummyCommand">
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

    is($dag->get_xml, $expected_xml, 'simple dag produces expected xml');
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
<operation name="top level" logDir="/tmp">
  <operationtype typeClass="Workflow::OperationType::Model">
    <inputproperty>some_external_input</inputproperty>
    <outputproperty>some_external_output</outputproperty>
  </operationtype>
  <operation name="some op">
    <operationtype typeClass="Workflow::OperationType::Command" lsfQueue="apipe" lsfResource="-M 25000000 -R 'select[mem&gt;25000] rusage[mem=25000]'" commandClass="Genome::WorkflowBuilder::Test::DummyCommand">
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
    is($dag->get_xml, $xml, 'xml round trip');
};

done_testing();
