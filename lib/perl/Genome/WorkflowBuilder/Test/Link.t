use strict;
use warnings;

use above 'Genome';
use Test::Exception;
use Test::More;


use_ok('Genome::WorkflowBuilder::Link');

subtest 'Typical Link' => sub {
    my $source_op = Genome::WorkflowBuilder::Command->create(
        name => 'source op',
        command => 'Genome::WorkflowBuilder::Test::DummyCommand'
    );
    my $destination_op = Genome::WorkflowBuilder::Command->create(
        name => 'destination op',
        command => 'Genome::WorkflowBuilder::Test::DummyCommand'
    );

    my $link = Genome::WorkflowBuilder::Link->create(
        source => $source_op, source_property => 'single_output',
        destination => $destination_op, destination_property => 'input'
    );

    my $expected_xml = '<link fromOperation="source op" fromProperty="single_output" toOperation="destination op" toProperty="input"/>';
    is($link->get_xml, $expected_xml, 'typical link produces expected xml');
};

subtest 'Parallel-By Link' => sub {
    my $source_op = Genome::WorkflowBuilder::Command->create(
        name => 'source op',
        command => 'Genome::WorkflowBuilder::Test::DummyCommand'
    );
    my $destination_op = Genome::WorkflowBuilder::Command->create(
        name => 'destination op',
        command => 'Genome::WorkflowBuilder::Test::DummyCommand',
        parallel_by => 'input',
    );

    my $link = Genome::WorkflowBuilder::Link->create(
        source => $source_op, source_property => 'many_output',
        destination => $destination_op, destination_property => 'input'
    );
    my $expected_xml = '<link fromOperation="source op" fromProperty="many_output" toOperation="destination op" toProperty="input"/>';
    is($link->get_xml, $expected_xml, 'parallelBy link produces expected xml');
};

subtest 'Input Connector' => sub {
    my $destination_op = Genome::WorkflowBuilder::Command->create(
        name => 'destination op',
        command => 'Genome::WorkflowBuilder::Test::DummyCommand'
    );

    my $link = Genome::WorkflowBuilder::Link->create(
        source_property => 'some_external_input',
        destination => $destination_op, destination_property => 'input'
    );

    my $expected_xml = '<link fromOperation="input connector" fromProperty="some_external_input" toOperation="destination op" toProperty="input"/>';
    is($link->get_xml, $expected_xml,
        'missing source operation uses input connector');
};

subtest 'Output Connector' => sub {
    my $source_op = Genome::WorkflowBuilder::Command->create(
        name => 'source op',
        command => 'Genome::WorkflowBuilder::Test::DummyCommand'
    );

    my $link = Genome::WorkflowBuilder::Link->create(
        source => $source_op, source_property => 'single_output',
        destination_property => 'some_external_output'
    );

    my $expected_xml = '<link fromOperation="source op" fromProperty="single_output" toOperation="output connector" toProperty="some_external_output"/>';
    is($link->get_xml, $expected_xml,
        'missing destination operation uses output connector');
};


subtest 'Valid Operation Type for Notify' => sub {
    dies_ok { Genome::WorkflowBuilder::Link->create(
        source => 'INVALID_OPERATION', source_property => 'foo',
        destination_property => 'bar') }, 'invalid source operation dies';

    dies_ok { Genome::WorkflowBuilder::Link->create(
        source_property => 'foo',
        destination => 'INVALID_OPERATION', destination_property => 'bar') },
        'invalid destination operation dies';
};

subtest 'Source Property Valid' => sub {
    my $source_op = Genome::WorkflowBuilder::Command->create(
        name => 'source op',
        command => 'Genome::WorkflowBuilder::Test::DummyCommand'
    );

    my $link = Genome::WorkflowBuilder::Link->create(
        source => $source_op, source_property => 'invalid_property',
        destination_property => 'some_external_output');

    eval {
        diag "Expect error message about invalid source property here:";
        $link->validate;
    };
    ok($@, 'invalid source property fails to validate');
};

subtest 'Destination Property Valid' => sub {
    my $destination_op = Genome::WorkflowBuilder::Command->create(
        name => 'destination op',
        command => 'Genome::WorkflowBuilder::Test::DummyCommand'
    );

    my $link = Genome::WorkflowBuilder::Link->create(
        source_property => 'some_external_input',
        destination => $destination_op,
        destination_property => 'invalid_property');

    eval {
        diag "Expect error message about invalid destination property here:";
        $link->validate;
    };
    ok($@, 'invalid destination property fails to validate');
};

subtest 'Parallel-By Data Source is_many' => sub {
    my $source_op = Genome::WorkflowBuilder::Command->create(
        name => 'source op',
        command => 'Genome::WorkflowBuilder::Test::DummyCommand'
    );
    my $destination_op = Genome::WorkflowBuilder::Command->create(
        name => 'destination op',
        command => 'Genome::WorkflowBuilder::Test::DummyCommand',
        parallel_by => 'input',
    );

    my $link = Genome::WorkflowBuilder::Link->create(
        source => $source_op, source_property => 'single_output',
        destination => $destination_op, destination_property => 'input'
    );

    eval {
        diag "Expect error message for destination property not being is_many";
        $link->validate;
    };

    ok($@, 'non is_many properties cannot be linked to parallelBy input');
};

# This tests that one can have a parallel_by operation with one singular input and one parallel_by is_many input
subtest 'Parallel-By Link Plus Single Link' => sub {
    my $source_op = Genome::WorkflowBuilder::Command->create(
        name => 'source op',
        command => 'Genome::WorkflowBuilder::Test::DummyCommand'
    );
    my $destination_op = Genome::WorkflowBuilder::Command->create(
        name => 'destination op',
        command => 'Genome::WorkflowBuilder::Test::DummyCommandTwoInputs',
        parallel_by => 'input',
    );

    my $link_many = Genome::WorkflowBuilder::Link->create(
        source => $source_op, source_property => 'many_output',
        destination => $destination_op, destination_property => 'input'
    );
    ok($link_many->validate, "parallelBy link validates");

    my $link_single = Genome::WorkflowBuilder::Link->create(
        source => $source_op, source_property => 'single_output',
        destination => $destination_op, destination_property => 'input_two'
    );
    ok($link_single->validate, "single link validates");

    my $expected_xml_many = '<link fromOperation="source op" fromProperty="many_output" toOperation="destination op" toProperty="input"/>';
    my $expected_xml_single = '<link fromOperation="source op" fromProperty="single_output" toOperation="destination op" toProperty="input_two"/>';

    is($link_many->get_xml, $expected_xml_many, 'parallelBy link produces expected xml');
    is($link_single->get_xml, $expected_xml_single, 'single link produces expected xml');
};

done_testing();
