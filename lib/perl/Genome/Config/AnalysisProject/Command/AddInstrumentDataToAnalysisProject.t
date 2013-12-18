#!/usr/bin/env genome-perl

use strict;
use warnings;

BEGIN {
    $ENV{UR_DBI_NO_COMMIT} = 1;
}

use above 'Genome';
use Genome::Test::Factory::InstrumentData::Solexa;

use Test::More;

my $class = 'Genome::Config::AnalysisProject::Command::AddInstrumentDataToAnalysisProject';
use_ok($class);

my $inst_data_1 = Genome::InstrumentData::Imported->create();
my $inst_data_2 = Genome::InstrumentData::Imported->create();

my $ap = Genome::Config::AnalysisProject->create(
    name => 'Test Project'
);

my $cmd = $class->create(
    instrument_data => [$inst_data_1, $inst_data_2],
    analysis_project => $ap
);

$cmd->execute();

for ($inst_data_1, $inst_data_2) {
    ok(Genome::Config::AnalysisProject::InstrumentDataBridge->get(
            analysis_project => $ap,
            instrument_data => $_,
        ), 'it should create a bridge entity for each instrument data');
}

done_testing();
