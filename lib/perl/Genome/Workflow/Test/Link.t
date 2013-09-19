use strict;
use warnings;

use above 'Genome';
use Test::More;


use_ok('Genome::Workflow::Link');

subtest 'Typical Link' => sub {
    my $source_op = Genome::Workflow::Command->create(
        name => 'source op',
        command => 'Genome::Workflow::Test::DummyCommand'
    );
    my $destination_op = Genome::Workflow::Command->create(
        name => 'destination op',
        command => 'Genome::Workflow::Test::DummyCommand'
    );

    my $link = Genome::Workflow::Link->create(
        source => $source_op, source_property => 'single_output',
        destination => $destination_op, destination_property => 'input'
    );

    my $expected_xml = '<link fromOperation="source op" fromProperty="single_output" toOperation="destination op" toProperty="input"/>';
    is($link->get_xml, $expected_xml, 'typical link produces expected xml');
};

subtest 'Parallel-By Link' => sub {
    my $source_op = Genome::Workflow::Command->create(
        name => 'source op',
        command => 'Genome::Workflow::Test::DummyCommand'
    );
    my $destination_op = Genome::Workflow::Command->create(
        name => 'destination op',
        command => 'Genome::Workflow::Test::DummyCommand',
        parallel_by => 'input',
    );

    my $link = Genome::Workflow::Link->create(
        source => $source_op, source_property => 'many_output',
        destination => $destination_op, destination_property => 'input'
    );
    my $expected_xml = '<link fromOperation="source op" fromProperty="many_output" toOperation="destination op" toProperty="input"/>';
    is($link->get_xml, $expected_xml, 'parallelBy link produces expected xml');
};

subtest 'Input Connector' => sub {
    my $destination_op = Genome::Workflow::Command->create(
        name => 'destination op',
        command => 'Genome::Workflow::Test::DummyCommand'
    );

    my $link = Genome::Workflow::Link->create(
        source_property => 'some_external_input',
        destination => $destination_op, destination_property => 'input'
    );

    my $expected_xml = '<link fromOperation="input connector" fromProperty="some_external_input" toOperation="destination op" toProperty="input"/>';
    is($link->get_xml, $expected_xml,
        'missing source operation uses input connector');
};

subtest 'Output Connector' => sub {
    my $source_op = Genome::Workflow::Command->create(
        name => 'source op',
        command => 'Genome::Workflow::Test::DummyCommand'
    );

    my $link = Genome::Workflow::Link->create(
        source => $source_op, source_property => 'single_output',
        destination_property => 'some_external_output'
    );

    my $expected_xml = '<link fromOperation="source op" fromProperty="single_output" toOperation="output connector" toProperty="some_external_output"/>';
    is($link->get_xml, $expected_xml,
        'missing destination operation uses output connector');
};


subtest 'Valid Operation Type' => sub {
    my $link_with_invalid_source = Genome::Workflow::Link->create(
        source => 'INVALID_OPERATION', source_property => 'foo',
        destination_property => 'bar'
    );
    eval {
        diag "Expect one error message about invalid source here:";
        $link_with_invalid_source->validate;
    };
    ok($@, 'invalid source operation fails to validate');

    my $link_with_invalid_destination = Genome::Workflow::Link->create(
        source_property => 'foo',
        destination => 'INVALID_OPERATION', destination_property => 'bar'
    );
    eval {
        diag "Expect one error message about invalid destination here:";
        $link_with_invalid_destination->validate;
    };
    ok($@, 'invalid destination operation fails to validate');
};

subtest 'Source Property Valid' => sub {
    my $source_op = Genome::Workflow::Command->create(
        name => 'source op',
        command => 'Genome::Workflow::Test::DummyCommand'
    );

    my $link = Genome::Workflow::Link->create(
        source => $source_op, source_property => 'invalid_property',
        destination_property => 'some_external_output');

    eval {
        diag "Expect error message about invalid source property here:";
        $link->validate;
    };
    ok($@, 'invalid source property fails to validate');
};

subtest 'Destination Property Valid' => sub {
    my $destination_op = Genome::Workflow::Command->create(
        name => 'destination op',
        command => 'Genome::Workflow::Test::DummyCommand'
    );

    my $link = Genome::Workflow::Link->create(
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
    my $source_op = Genome::Workflow::Command->create(
        name => 'source op',
        command => 'Genome::Workflow::Test::DummyCommand'
    );
    my $destination_op = Genome::Workflow::Command->create(
        name => 'destination op',
        command => 'Genome::Workflow::Test::DummyCommand',
        parallel_by => 'input',
    );

    my $link = Genome::Workflow::Link->create(
        source => $source_op, source_property => 'single_output',
        destination => $destination_op, destination_property => 'input'
    );

    eval {
        diag "Expect error message for destination property not being is_many";
        $link->validate;
    };

    ok($@, 'non is_many properties cannot be linked to parallelBy input');
};

done_testing();
