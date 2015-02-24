#!/usr/bin/env genome-perl

use strict;
use warnings;

BEGIN {
    $ENV{UR_DBI_NO_COMMIT} = 1;
};

use Test::More tests => 5;

use above 'Genome';

use Genome::Test::Factory::InstrumentData::Solexa;

my $command_class = 'Genome::Config::AnalysisProject::Command::Reprocess';
use_ok($command_class);

my $analysis_project = Genome::Config::AnalysisProject->create(
    name => 'Test Project for Reprocess Command',
);

my $bridge_class = 'Genome::Config::AnalysisProject::InstrumentDataBridge';
my $possible_statuses = $bridge_class->__meta__->property(property_name => 'status')->valid_values;

my @instrument_data;
my @bridges;

for my $i (0..$#$possible_statuses) {
    my $instrument_data = Genome::Test::Factory::InstrumentData::Solexa->setup_object(lane => $i);
    push @instrument_data, $instrument_data;

    my $bridge = $bridge_class->create(
        analysis_project => $analysis_project,
        instrument_data => $instrument_data,
        status => $possible_statuses->[$i],
    );
    push @bridges, $bridge;
}

my %statuses = map(($_->status => 1), @bridges);
cmp_ok(scalar(keys %statuses), '>', 1, 'bridges do not all have the same initial status');

my $cmd = $command_class->create(analysis_project => $analysis_project);
isa_ok($cmd, $command_class, 'created command');

ok($cmd->execute, 'executed command');

subtest 'all bridges are rescheduled' => sub {
    for my $bridge (@bridges) {
        is($bridge->status, 'new', 'bridge is rescheduled');
    }
};
