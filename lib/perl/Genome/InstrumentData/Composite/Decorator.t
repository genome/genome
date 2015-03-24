#!/usr/bin/env genome-perl

use strict;
use warnings;

BEGIN {
    $ENV{UR_DBI_NO_COMMIT} = 1;
};

use Test::More tests => 3;
use above 'Genome';
use Workflow;

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
        my($class, $operation, $params) = @_;
        $call_count++;
        $called_params = $params;
    }
}

my $test_params = 'test params';
my $decoration = {
    name => 'test-decorator',
    params => $test_params,
};

my $operation = Workflow::Operation->__define__(name => 'test operation for decorators');

Genome::InstrumentData::Composite::Decorator->decorate($operation, $decoration);

is($call_count, 1, 'decorator got called');
is($called_params, $test_params, 'decorator was passed params');
