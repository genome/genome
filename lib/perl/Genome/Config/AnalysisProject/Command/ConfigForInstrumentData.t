#!/usr/bin/env genome-perl

use strict;
use warnings;

BEGIN {
    $ENV{UR_DBI_NO_COMMIT} = 1;
    $ENV{UR_COMMAND_DUMP_STATUS_MESSAGES} = 1;
    $ENV{UR_COMMAND_DUMP_DEBUG_MESSAGES} = 1;
};

use Test::More tests => 7;

use above 'Genome';

use Genome::Test::Factory::InstrumentData::Solexa;

my $class = 'Genome::Config::AnalysisProject::Command::ConfigForInstrumentData';
use_ok($class) or die;

my $ap = Genome::Config::AnalysisProject->__define__(
    name => '__TEST_AnP__',
    status => 'In Progress',
);
ok($ap, 'define AnP');

my $instrument_data = Genome::Test::Factory::InstrumentData::Solexa->setup_object;
ok($instrument_data, 'define instdata');

#test with missing Analysis Project
my $command = $class->execute( instrument_data => [$instrument_data] );
ok($command->result, 'execute w/o ANp instdata bridge');
like($command->warning_message, qr/No analysis-projects associated with instrument data/, 'correct warning_message');

#connect AnP to instdata
my $bridge = Genome::Config::AnalysisProject::InstrumentDataBridge->__define__(
    analysis_project => $ap,
    instrument_data => $instrument_data,
);
ok($bridge, 'define AnP instdata bridge');

$command = $class->execute(instrument_data => [$instrument_data]);
ok($command->result, 'execute');

done_testing();
