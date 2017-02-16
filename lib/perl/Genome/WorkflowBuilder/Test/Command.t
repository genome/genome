use strict;
use warnings;

use above 'Genome';
use Test::More;

use Genome::Test::Config qw(setup_config);


use_ok('Genome::WorkflowBuilder::Command');

my $lsf_queue_build_worker_alt = Genome::Config::get('lsf_queue_build_worker_alt');

subtest 'Typical Command' => sub {
    my $op = Genome::WorkflowBuilder::Command->create(
        name => 'some op',
        command => 'Genome::WorkflowBuilder::Test::DummyCommand'
    );

    my $expected_xml = <<EOS;
<?xml version="1.0"?>
<operation name="some op">
  <operationtype typeClass="Workflow::OperationType::Command" commandClass="Genome::WorkflowBuilder::Test::DummyCommand" lsfQueue="$lsf_queue_build_worker_alt" lsfResource="-M 25000000 -R 'select[mem&gt;25000] rusage[mem=25000]'">
    <inputproperty>input</inputproperty>
    <outputproperty>many_output</outputproperty>
    <outputproperty>result</outputproperty>
    <outputproperty>single_output</outputproperty>
  </operationtype>
</operation>
EOS

    is($op->get_xml, $expected_xml, 'typical command produces expected xml');
};

subtest 'Command with different attributes in xml than on Command class' => sub {

    my $expected_xml = <<EOS;
<?xml version="1.0"?>
<operation name="some op">
  <operationtype typeClass="Workflow::OperationType::Command" commandClass="Genome::WorkflowBuilder::Test::DummyCommand" lsfQueue="bob" lsfResource="-M 25000 -R 'select[mem&gt;25] rusage[mem=25]'">
    <inputproperty>input</inputproperty>
    <outputproperty>many_output</outputproperty>
    <outputproperty>result</outputproperty>
    <outputproperty>single_output</outputproperty>
  </operationtype>
</operation>
EOS

    my $op = Genome::WorkflowBuilder::Command->from_xml($expected_xml);
    is($op->get_xml, $expected_xml, 'Xml attributes take priority');
};

subtest 'Parallel-By Command' => sub {
    my $op = Genome::WorkflowBuilder::Command->create(
        name => 'some op',
        command => 'Genome::WorkflowBuilder::Test::DummyCommand',
        parallel_by => 'input',
    );

    my $expected_xml = <<EOS;
<?xml version="1.0"?>
<operation name="some op" parallelBy="input">
  <operationtype typeClass="Workflow::OperationType::Command" commandClass="Genome::WorkflowBuilder::Test::DummyCommand" lsfQueue="$lsf_queue_build_worker_alt" lsfResource="-M 25000000 -R 'select[mem&gt;25000] rusage[mem=25000]'">
    <inputproperty>input</inputproperty>
    <outputproperty>many_output</outputproperty>
    <outputproperty>result</outputproperty>
    <outputproperty>single_output</outputproperty>
  </operationtype>
</operation>
EOS

    is($op->get_xml, $expected_xml, 'parallelBy command produces expected xml');
};

subtest 'Invalid Parallel-By Command' => sub {
    my $op = Genome::WorkflowBuilder::Command->create(
        name => 'some op',
        command => 'Genome::WorkflowBuilder::Test::DummyCommand',
        parallel_by => 'invalid_input',
    );

    eval {
        diag "Expect one error message about failing to validate inputs here:";
        $op->validate;
    };

    ok($@, "parallelBy on non input_property doesn't validate");
};

subtest 'Invalid Command Name' => sub {
    my $op = Genome::WorkflowBuilder::Command->create(
        name => 'input connector',
        command => 'Genome::WorkflowBuilder::Test::DummyCommand',
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
<operation name="some op" parallelBy="input">
  <operationtype typeClass="Workflow::OperationType::Command" commandClass="Genome::WorkflowBuilder::Test::DummyCommand" lsfQueue="$lsf_queue_build_worker_alt" lsfResource="-M 25000000 -R 'select[mem&gt;25000] rusage[mem=25000]'">
    <inputproperty>input</inputproperty>
    <outputproperty>many_output</outputproperty>
    <outputproperty>result</outputproperty>
    <outputproperty>single_output</outputproperty>
  </operationtype>
</operation>
EOS

    my $op = Genome::WorkflowBuilder::Command->from_xml($xml);
    is($op->get_xml, $xml, 'xml round trip');
};

subtest 'unspecified_operation_type_attributes' => sub {
    my $op = Genome::WorkflowBuilder::Command->create(
        name => 'some op',
        command => 'Genome::WorkflowBuilder::Test::DummyCommand',
    );
    my %got = $op->operation_type_attributes;
    my %expected = (
        lsfQueue => $lsf_queue_build_worker_alt,
        lsfResource => "-M 25000000 -R 'select[mem>25000] rusage[mem=25000]'",
        commandClass => 'Genome::WorkflowBuilder::Test::DummyCommand',
    );
    is_deeply(\%got, \%expected, 'got attributes from command');
};

subtest 'specified_operation_type_attributes' => sub {
    my $op = Genome::WorkflowBuilder::Command->create(
        name => 'some op',
        command => 'Genome::WorkflowBuilder::Test::DummyCommand',
    );
    my $queue = 'not apipe';
    $op->lsf_queue($queue);
    $op->lsf_project('specified project');
    my %got = $op->operation_type_attributes;
    my %expected = (
        lsfQueue => $queue,
        lsfResource => "-M 25000000 -R 'select[mem>25000] rusage[mem=25000]'",
        lsfProject => 'specified project',
        commandClass => 'Genome::WorkflowBuilder::Test::DummyCommand',
    );
    is_deeply(\%got, \%expected, 'got attributes as specified first, then from command');
};

subtest 'command with config (calculated_default)' => sub {
    my ($temp_dirs, $new_temp_dir) = Genome::Test::Config::temp_dir_helper();
    local $ENV{XGENOME_CONFIG_SNAP} = $new_temp_dir->();
    local $ENV{XGENOME_CONFIG_HOME} = $new_temp_dir->();
    local $ENV{XGENOME_CONFIG_DIRS} = $new_temp_dir->();

    setup_config(
        spec => {
            dummy_lsf_resource => {},
            dummy_lsf_queue => {},
        },
        global => {
            dummy_lsf_resource => 'my_dummy_lsf_resource',
            dummy_lsf_queue => 'my_dummy_lsf_queue',
        },
    );
    local $INC{'Genome/WorkflowBuilder/Test/DummyCommandCalculatedDefault.pm'} = 1;
    UR::Object::Type->define(
        class_name => 'Genome::WorkflowBuilder::Test::DummyCommandCalculatedDefault',
        is => [qw(Command Genome::Configurable)],

        has_input => ['input'],

        has_output => [
            many_output => {
                is_many => 1,
            },
            single_output => { },
        ],

        has_param => [
            lsf_resource => {
                config => 'dummy_lsf_resource',
            },
            lsf_queue => {
                config => 'dummy_lsf_queue',
            },
        ],
    );

    my $op = Genome::WorkflowBuilder::Command->create(
        name => 'some op',
        command => 'Genome::WorkflowBuilder::Test::DummyCommandCalculatedDefault'
    );

    my $expected_xml = <<EOS;
<?xml version="1.0"?>
<operation name="some op">
  <operationtype typeClass="Workflow::OperationType::Command" commandClass="Genome::WorkflowBuilder::Test::DummyCommandCalculatedDefault" lsfQueue="my_dummy_lsf_queue" lsfResource="my_dummy_lsf_resource">
    <inputproperty>input</inputproperty>
    <outputproperty>many_output</outputproperty>
    <outputproperty>result</outputproperty>
    <outputproperty>single_output</outputproperty>
  </operationtype>
</operation>
EOS

    is($op->get_xml, $expected_xml, 'typical command produces expected xml');
};

done_testing();
