#!/usr/bin/env genome-perl

use strict;
use warnings;

BEGIN {
    $ENV{UR_DBI_NO_COMMIT} = 1;
};

use Test::More tests => 3;
use above 'Genome';

my $pkg = 'Genome::InstrumentData::Composite::Decorator';
use_ok($pkg);

my $call_count = 0;
my $called_params;
{
    package Genome::InstrumentData::Composite::Decorator::TestDecorator;
    class Genome::InstrumentData::Composite::Decorator::TestDecorator {
        is => 'Genome::InstrumentData::Composite::Decorator::Base',
    };
    sub decorate {
        my($class, $operation, $workflow, $params) = @_;
        $call_count++;
        $called_params = $params;
    }
}

my $test_params = 'test params';
my $decoration = {
    name => 'test-decorator',
    params => $test_params,
};

my $dag = Genome::WorkflowBuilder::DAG->create(name => 'test workflow');
my $operation = Genome::WorkflowBuilder::Command->create(name => 'test operation for decorators', command => 'Genome::Model::Tools::Example1');
$dag->add_operation($operation);

Genome::InstrumentData::Composite::Decorator->decorate($operation, $dag, $decoration);

is($call_count, 1, 'decorator got called');
is($called_params, $test_params, 'decorator was passed params');
