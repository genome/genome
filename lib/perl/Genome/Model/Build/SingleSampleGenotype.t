#!/usr/bin/env genome-perl

use strict;
use warnings;

BEGIN {
    $ENV{UR_DBI_NO_COMMIT} = 1;
};

use above 'Genome';

use Genome::Test::Factory::InstrumentData::Solexa;
use Genome::Test::Factory::Model::SingleSampleGenotype;

use Test::More tests => 3;

my $pkg = 'Genome::Model::Build::SingleSampleGenotype';
use_ok($pkg);

my $instrument_data = Genome::Test::Factory::InstrumentData::Solexa->setup_object();
my $model = Genome::Test::Factory::Model::SingleSampleGenotype->setup_object();
$model->add_instrument_data($instrument_data);

my $build = $pkg->create(
    model_id => $model->id,
);
isa_ok($build, $pkg, 'created build');

is($build->instrument_data, $instrument_data, 'build has instrument data assigned');
