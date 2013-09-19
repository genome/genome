use strict;
use warnings;

use above 'Genome';
use Test::More;


use_ok('Genome::Workflow::Command');

subtest 'Typical Command' => sub {
    my $op = Genome::Workflow::Command->create(
        name => 'some op',
        command => 'Genome::Workflow::Test::DummyCommand'
    );

    my $expected_xml = <<EOS;
<?xml version="1.0"?>
<operation name="some op"><operationtype typeClass="Workflow::OperationType::Command" lsfQueue="apipe" lsfResource="-M 25000000 -R 'select[mem&gt;25000] rusage[mem=25000]'" commandClass="Genome::Workflow::Test::DummyCommand"><inputproperty>input</inputproperty><outputproperty>many_output</outputproperty><outputproperty>result</outputproperty><outputproperty>single_output</outputproperty></operationtype></operation>
EOS

    is($op->get_xml, $expected_xml, 'typical command produces expected xml');
};

subtest 'Parallel-By Command' => sub {
    my $op = Genome::Workflow::Command->create(
        name => 'some op',
        command => 'Genome::Workflow::Test::DummyCommand',
        parallel_by => 'input',
    );

    my $expected_xml = <<EOS;
<?xml version="1.0"?>
<operation name="some op" parallelBy="input"><operationtype typeClass="Workflow::OperationType::Command" lsfQueue="apipe" lsfResource="-M 25000000 -R 'select[mem&gt;25000] rusage[mem=25000]'" commandClass="Genome::Workflow::Test::DummyCommand"><inputproperty>input</inputproperty><outputproperty>many_output</outputproperty><outputproperty>result</outputproperty><outputproperty>single_output</outputproperty></operationtype></operation>
EOS

    is($op->get_xml, $expected_xml, 'parallelBy command produces expected xml');
};

subtest 'Invalid Parallel-By Command' => sub {
    my $op = Genome::Workflow::Command->create(
        name => 'some op',
        command => 'Genome::Workflow::Test::DummyCommand',
        parallel_by => 'invalid_input',
    );

    eval {
        diag "Expect one error message about failing to validate inputs here:";
        $op->validate;
    };

    ok($@, "parallelBy on non input_property doesn't validate");
};

subtest 'Invalid Command Name' => sub {
    my $op = Genome::Workflow::Command->create(
        name => 'input connector',
        command => 'Genome::Workflow::Test::DummyCommand',
    );

    eval {
        diag "Expect one error message about invalid operation name:";
        $op->validate;
    };

    ok($@, 'invalid operation name fails to validate');
};

subtest 'XML Round Trip' => sub {
    my $xml = <<EOS;
<?xml version="1.0"?>
<operation name="some op"><operationtype typeClass="Workflow::OperationType::Command" lsfQueue="apipe" lsfResource="-M 25000000 -R 'select[mem&gt;25000] rusage[mem=25000]'" commandClass="Genome::Workflow::Test::DummyCommand"><inputproperty>input</inputproperty><outputproperty>many_output</outputproperty><outputproperty>result</outputproperty><outputproperty>single_output</outputproperty></operationtype></operation>
EOS

    my $op = Genome::Workflow::Command->from_xml($xml);
    is($op->get_xml, $xml, 'xml round trip');
};

done_testing();
