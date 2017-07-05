#!/usr/bin/env genome-perl

use strict;
use warnings;

BEGIN {
    $ENV{UR_DBI_NO_COMMIT} = 1;
}

use above 'Genome';
use Genome::Test::Factory::InstrumentData::Solexa;

use Test::More tests => 6;

my $class = 'Genome::Model::CwlPipeline::Command::Input::Remove';
use_ok($class);

my $value = Genome::Test::Factory::InstrumentData::Solexa->setup_object;
my $other_value = Genome::Test::Factory::InstrumentData::Solexa->setup_object;

my $m = Genome::Model::CwlPipeline->__define__(
    name => 'cwl pipeline input remove test model',
);

for my $v ($value, $other_value) {
    my $input = Genome::Model::Input->create(
        name => 'instrument_data',
        value_id => $v->id,
        value_class_name => $v->class,
        model_id => $m->id,
    );
}

my @prior_inputs = $m->inputs;
is(scalar(@prior_inputs), 2, 'inputs initially exist');

my $remove_cmd = $class->create(
    model => $m,
    name => 'instrument_data',
    value => $value->id,
);
isa_ok($remove_cmd, $class, 'created command');

ok($remove_cmd->execute, 'executed command');

my @inputs = $m->inputs;
is(scalar(@inputs), 1, 'removed one input');
is($inputs[0]->value, $other_value, 'other input remains');
