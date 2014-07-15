#!/usr/bin/env genome-perl

use strict;
use warnings;

BEGIN {
    $ENV{UR_DBI_NO_COMMIT} = 1;
};

use Test::More tests => 3;

use above 'Genome';

use Genome::Test::Factory::InstrumentData::Solexa;

my $class = 'Genome::Config::AnalysisProject::Command::ConfigForInstrumentData';

use_ok($class);

my $instrument_data = Genome::Test::Factory::InstrumentData::Solexa->setup_object;

#test with missing Analysis Project
my $command = $class->create( instrument_data => [$instrument_data] );
isa_ok($command, $class, 'created command');
ok($command->execute, 'executed command');

