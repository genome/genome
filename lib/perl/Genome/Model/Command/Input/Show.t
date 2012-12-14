#!/usr/bin/env genome-perl

use strict;
use warnings;

use above 'Genome';

use Data::Dumper 'Dumper';
use Test::More;

use_ok('Genome::Model::Command::Input::Show') or die;

my $show = Genome::Model::Command::Input::Show->create();
ok($show, 'create');
ok(!$show->execute, 'execute');

done_testing();
