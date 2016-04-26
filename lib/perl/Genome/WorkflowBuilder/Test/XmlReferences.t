#!/usr/bin/env genome-perl

use strict;
use warnings;

BEGIN {
    $ENV{UR_DBI_NO_COMMIT} = 1;
}

use above 'Genome';

use Sub::Install qw();

use Test::More tests => 6;

my $test_dir = __FILE__ . '.d';

use Genome::WorkflowBuilder::Test::DummyCommand;

my $count = 0;
my $expected = 'expected output';
Sub::Install::install_sub({
    code => sub {
        my $self = shift;
        $self->single_output($expected);

        return ++$count;
    },
    into => 'Genome::WorkflowBuilder::Test::DummyCommand',
    as   => '_execute_body',
});

#Since the outputs are "required" properties, set defaults here
my $meta = Genome::WorkflowBuilder::Test::DummyCommand->__meta__;
$meta->property(property_name => 'single_output')->default_value('unexpected output');
$meta->property(property_name => 'many_output')->default_value([qw(not used)]);

my $xml = File::Spec->join($test_dir, 'outer.xml');
my $dag = Genome::WorkflowBuilder::DAG->from_xml_filename($xml);
isa_ok($dag, 'Genome::WorkflowBuilder::DAG', 'constructed DAG from XML');

my $inner = $dag->operation_named('inner');
isa_ok($inner, 'Genome::WorkflowBuilder::DAG', 'created inner DAG with name from outer DAG');

my $dummy = $inner->operation_named('dummy');
is($dummy->command, 'Genome::WorkflowBuilder::Test::DummyCommand', 'created operation');

my $outputs = $dag->execute_inline({input => 'test'});
ok($outputs, 'executed DAG');
is($outputs->{single_output}, $expected, 'got output expected');
is($count, 1, 'command was executed as expected');

done_testing();
