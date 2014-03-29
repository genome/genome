#!/usr/bin/env genome-perl

use strict;
use warnings;

use above 'Genome';

use Test::More;

use_ok('Genome::Model::GenotypeMicroarray::Command::Alleles') or die;
use_ok('Genome::Model::GenotypeMicroarray::Test') or die;

my $build = Genome::Model::GenotypeMicroarray::Test->example_build;
my $cmd = Genome::Model::GenotypeMicroarray::Command::Alleles->create(build => $build);
ok($cmd, 'create alleels command');
ok($cmd->execute, 'execute alleles command');
is_deeply($cmd->alleles, {AA => 1, AG => 3, CC => 2, GG => 2, TC => 1, TT => 1, }, 'alleles');
is($cmd->total_genotypes, 10, 'total genotypes');

done_testing();
