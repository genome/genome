#!/usr/bin/env genome-perl

use strict;
use warnings;

use above 'Genome';
use Test::More tests => 2;

use_ok('Genome::Model::Tools::Lims');

my $cmd = Genome::Model::Tools::Lims->create(args => ['--help']);
ok($cmd->execute(), '`gmt lims --help` works');
