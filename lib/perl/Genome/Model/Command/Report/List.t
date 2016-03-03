#!/usr/bin/env genome-perl

use strict;
use warnings;

use Test::More;

use_ok('Genome::Model::Command::Report::List') or die;

ok(Genome::Model::Command::Report::List->execute()->result, 'execute');
ok(Genome::Model::Command::Report::List->execute(type_names => 1)->result, 'execute type_names => 1');
ok(Genome::Model::Command::Report::List->execute(type_name => 'reference alignment')->result, 'execute type_name => reference alignment');
ok(Genome::Model::Command::Report::List->execute(generic => 1)->result, 'execute generic => 1');

my $list = Genome::Model::Command::Report::List->create(type_name => 'no way this is a model');
ok($list, 'create w/ invalid type name');
ok(!$list->execute, 'execute failed w/ invalid type name');
$list = Genome::Model::Command::Report::List->create(type_name => 'reference alignment', all => 1);
ok($list, 'create w/ more than one list method');
ok(!$list->execute, 'execute failed w/ more than one list method');

done_testing();
