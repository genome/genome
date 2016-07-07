#!/usr/bin/env genome-perl

use strict;
use warnings;


use Sub::Override;
use Test::More tests => 3;

use above 'Genome';
use Genome::Test::Factory::AnalysisProject;

my $class = 'Genome::SoftwareResult::Command::DiskUsage';

use_ok($class);

my $cmd = $class->create();
isa_ok($cmd, $class, 'created command instance with default parameters');

my $anp = Genome::Test::Factory::AnalysisProject->setup_object(name => 'Disk Usage Test AnP');

my $override = Sub::Override->new(
    "${class}::allocation_summary" => sub {
        return [[$anp->id, '100' * (1024 ** 3)], ['nobody@example.com', '50' * (1024 ** 3)]];
    },
);

ok($cmd->execute, 'executed command (skipping SQL query)');
