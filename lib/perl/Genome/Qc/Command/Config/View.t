#!/usr/bin/env genome-perl

use strict;
use warnings;

use above "Genome";

use Genome::Test::Factory::Qc::Config;

use Test::More tests => 3;

my $class = 'Genome::Qc::Command::Config::View';
use_ok($class) or die;

my $qc = Genome::Test::Factory::Qc::Config->setup_object(config => '{"key1": "value1"}', name => 'view_test', type => 'all');
my $cmd = $class->create(qc_config => $qc);
isa_ok($cmd, $class, 'created view command');
ok($cmd->execute, 'executed view command');

done_testing;
