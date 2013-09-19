use strict;
use warnings;

use above 'Genome';
use Test::More;


use_ok('Genome::Workflow::Link');

test_typical_link();

test_parallel_by_link();

test_input_connector();
test_output_connector();

test_validate_operations_are_valid_type();
test_validate_source_property();
test_validate_destination_property();
test_validate_parallel_by_source_is_many();

done_testing();


sub test_typical_link {
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

    my $expected_xml = '<link leftOperation="source op" leftProperty="single_output" rightOperation="destination op" rightProperty="input"/>';
    is($link->get_xml, $expected_xml, 'typical link produces expected xml');
}

sub test_parallel_by_link {
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
    my $expected_xml = '<link leftOperation="source op" leftProperty="many_output" rightOperation="destination op" rightProperty="input"/>';
    is($link->get_xml, $expected_xml, 'parallelBy link produces expected xml');
}

sub test_input_connector {
    my $destination_op = Genome::Workflow::Command->create(
        name => 'destination op',
        command => 'Genome::Workflow::Test::DummyCommand'
    );

    my $link = Genome::Workflow::Link->create(
        source_property => 'some_external_input',
        destination => $destination_op, destination_property => 'input'
    );

    my $expected_xml = '<link leftOperation="input connector" leftProperty="some_external_input" rightOperation="destination op" rightProperty="input"/>';
    is($link->get_xml, $expected_xml,
        'missing source operation uses input connector');
}

sub test_output_connector {
    my $source_op = Genome::Workflow::Command->create(
        name => 'source op',
        command => 'Genome::Workflow::Test::DummyCommand'
    );

    my $link = Genome::Workflow::Link->create(
        source => $source_op, source_property => 'single_output',
        destination_property => 'some_external_output'
    );

    my $expected_xml = '<link leftOperation="source op" leftProperty="single_output" rightOperation="output connector" rightProperty="some_external_output"/>';
    is($link->get_xml, $expected_xml,
        'missing destination operation uses output connector');
}


sub test_validate_operations_are_valid_type {
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
}

sub test_validate_source_property {
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
}

sub test_validate_destination_property {
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
}

sub test_validate_parallel_by_source_is_many {
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
}
