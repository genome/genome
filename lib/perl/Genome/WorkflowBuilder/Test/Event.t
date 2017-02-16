use strict;
use warnings;

use above 'Genome';
use Test::More;


use_ok('Genome::WorkflowBuilder::Event');

subtest 'Typical Event' => sub {
    my $op = Genome::WorkflowBuilder::Event->create(
        name => 'some event',
        event_id => '12345',
    );

    my $expected_xml = <<EOS;
<?xml version="1.0"?>
<operation name="some event">
  <operationtype typeClass="Workflow::OperationType::Event" eventId="12345"/>
</operation>
EOS

    is($op->get_xml, $expected_xml, 'typical command produces expected xml');
};

subtest 'Event roundtrip with attributes' => sub {
    my $expected_xml = <<EOS;
<?xml version="1.0"?>
<operation name="some other event">
  <operationtype typeClass="Workflow::OperationType::Event" eventId="12345" lsfProject="foo-bar" lsfQueue="bob" lsfResource="baz">
    <inputproperty>some_input</inputproperty>
    <outputproperty>some_output</outputproperty>
  </operationtype>
</operation>
EOS

    my $op = Genome::WorkflowBuilder::Event->from_xml($expected_xml);
    is($op->get_xml, $expected_xml, 'xml_roundtrip with attributes and inputs/outputs');
};

done_testing();
