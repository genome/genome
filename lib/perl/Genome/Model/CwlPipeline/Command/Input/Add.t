#!/usr/bin/env genome-perl

use strict;
use warnings;

BEGIN {
    $ENV{UR_DBI_NO_COMMIT} = 1;
}

use above 'Genome';
use Genome::Test::Factory::InstrumentData::Solexa;

use Test::More tests => 5;

my $class = 'Genome::Model::CwlPipeline::Command::Input::Add';
use_ok($class);

my $value = Genome::Test::Factory::InstrumentData::Solexa->setup_object;

my $m = Genome::Model::CwlPipeline->__define__(
    name => 'cwl pipeline input add test model',
);


my $add_cmd = $class->create(
    model => $m,
    name => 'instrument_data',
    value => $value->id,
);
isa_ok($add_cmd, $class, 'created command');

ok($add_cmd->execute, 'executed command');

my @inputs = $m->inputs;
is(scalar(@inputs), 1, 'added an input');
is($inputs[0]->value, $value, 'input has correct value');
