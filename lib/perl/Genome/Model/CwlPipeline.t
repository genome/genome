#!/usr/bin/env genome-perl

use strict;
use warnings;

BEGIN {
    $ENV{UR_DBI_NO_COMMIT} = 1;
}

use above 'Genome';
use Genome::Test::Factory::InstrumentData::Solexa;

use Test::More tests => 5;

my $class = 'Genome::Model::CwlPipeline';
use_ok($class);

my $instrument_data = Genome::Test::Factory::InstrumentData::Solexa->setup_object;

my $pp = Genome::ProcessingProfile::CwlPipeline->__define__(id => -100);

my $model = $class->create(
    processing_profile_id => $pp->id,
    subject => $instrument_data->sample,
    input_data => {
        instrument_data => [$instrument_data->id],
        awesomeness     => 'lots',
    },
);
isa_ok($model, $class, 'created model');

my @inputs = $model->inputs;
is(scalar(@inputs), 2, 'created inputs with model');

my ($ii) = grep { $_->name eq 'instrument_data' } @inputs;
is($ii->value, $instrument_data, 'set instrument data input to object');

my ($ai) = grep { $_->name eq 'awesomeness' } @inputs;
is($ai->value, 'lots', 'created text input');
