#!/usr/bin/env genome-perl

use strict;
use warnings;

use above "Genome";
use Genome::Test::Factory::Qc::Config;

use Test::More tests => 3;

my $class = 'Genome::Qc::Command::Config::Diff';
use_ok($class) or die;

my $qc1 = Genome::Test::Factory::Qc::Config->setup_object(config => '{"key": "value1"}', name => 'diff_test_a', type => 'all');
my $qc2 = Genome::Test::Factory::Qc::Config->setup_object(config => '{"key": "value2"}', name => 'diff_test_b', type => 'all');

my $cmd = $class->create(qc_config_a => $qc1, qc_config_b => $qc2);
isa_ok($cmd, $class, 'created diff command');
ok($cmd->execute, 'executed diff command');

done_testing;
