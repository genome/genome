#!/usr/bin/env genome-perl

BEGIN {
    $ENV{UR_DBI_NO_COMMIT} = 1;
    $ENV{UR_COMMAND_DUMP_STATUS_MESSAGES} = 1;
    $ENV{UR_COMMAND_DUMP_DEBUG_MESSAGES} = 1;
}

use strict;
use warnings;

use above 'Genome';

require Genome::Utility::Test;
use Test::More;

my $class = 'Genome::InstrumentData::Command::Import::WorkFlow::AddProcessToInstrumentData';
use_ok($class) or die;

my $process = Genome::InstrumentData::Command::Import::Process->__define__;
ok($process, 'define process');

my @instdata = map { Genome::InstrumentData::Imported->__define__(id => $_); } (qw/ 111 112 113 /);
is(@instdata, 3, 'define instdata');

my $cmd = $class->execute(
    instrument_data => \@instdata,
    process => $process,
);
ok($cmd->result, 'execute');

my @instdata_from_process = $process->instrument_data;
is_deeply(\@instdata_from_process, \@instdata, 'got correct instdata');

done_testing();
