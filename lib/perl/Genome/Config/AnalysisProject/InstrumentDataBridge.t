#!/usr/bin/env genome-perl

use strict;
use warnings;

use above 'Genome';

use Genome::Test::Factory::AnalysisProject;
use Genome::Test::Factory::InstrumentData::Solexa;
use Test::More;

use_ok('Genome::Config::AnalysisProject::InstrumentDataBridge') or die;

my $bridge = Genome::Config::AnalysisProject::InstrumentDataBridge->create(
    analysis_project => Genome::Test::Factory::AnalysisProject->setup_object,
    instrument_data => Genome::Test::Factory::InstrumentData::Solexa->setup_object,
    status => 'New',
);
ok($bridge, 'create bridge') or die;
is($bridge->__display_name__, sprintf('bridge between instrument data (%s) and analysis project (%s)', $bridge->instrument_data->id, $bridge->analysis_project->id), '__display_name__');

done_testing();
