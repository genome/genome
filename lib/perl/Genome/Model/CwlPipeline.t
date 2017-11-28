#!/usr/bin/env genome-perl

use strict;
use warnings;

BEGIN {
    $ENV{UR_DBI_NO_COMMIT} = 1;
}

use above 'Genome';
use Genome::Test::Factory::InstrumentData::Solexa;
use Genome::Test::Factory::Individual;

use Test::More tests => 7;

my $class = 'Genome::Model::CwlPipeline';
use_ok($class);

my $instrument_data = Genome::Test::Factory::InstrumentData::Solexa->setup_object;
my $arbitrary_object = Genome::Test::Factory::Individual->setup_object;

my $pp = Genome::ProcessingProfile::CwlPipeline->__define__(id => -100);

my @params = (
    processing_profile_id => $pp->id,
    subject => $instrument_data->sample,
    input_data => {
        instrument_data => [$instrument_data->id],
        awesomeness     => 'lots',
        arbitrary_object => {
            value_class_name => $arbitrary_object->class,
            value_id => $arbitrary_object->id,
        },
    },
);
my $model = $class->create(@params);
isa_ok($model, $class, 'created model');

my @inputs = $model->inputs;
is(scalar(@inputs), 3, 'created inputs with model');

my ($ii) = grep { $_->name eq 'instrument_data' } @inputs;
is($ii->value, $instrument_data, 'set instrument data input to object');

my ($ai) = grep { $_->name eq 'awesomeness' } @inputs;
is($ai->value, 'lots', 'created text input');

my ($ao) = grep { $_->name eq 'arbitrary_object' } @inputs;
is($ao->value, $arbitrary_object, 'set arbitrary object as in put');

my $model2 = $class->get(@params);
is($model2, $model, 'got our model that we created');
