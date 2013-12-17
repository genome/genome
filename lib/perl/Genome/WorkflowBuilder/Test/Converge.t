use strict;
use warnings;

use above 'Genome';
use Test::More;


subtest 'XML Round Trip' => sub {
    my $xml = <<EOS;
<?xml version="1.0"?>
<operation name="some op">
  <operationtype typeClass="Workflow::OperationType::Converge">
    <inputproperty>input_01</inputproperty>
    <inputproperty>input_02</inputproperty>
    <inputproperty>input_03</inputproperty>
    <outputproperty>converge_output</outputproperty>
    <outputproperty>result</outputproperty>
  </operationtype>
</operation>
EOS

    my $op = Genome::WorkflowBuilder::Converge->from_xml($xml);
    is($op->get_xml, $xml, 'xml round trip');
};


done_testing();
