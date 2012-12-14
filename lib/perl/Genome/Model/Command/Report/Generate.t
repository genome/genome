#!/usr/bin/env genome-perl

use strict;
use warnings;

use above 'Genome';

require File::Temp;
use Test::More;

use_ok('Genome::Model::Command::Report::Generate') or die;

my $build = Genome::Model::Build->get(107664200); # build for apipe-test-03-MC16s
ok($build, 'Got MC16s build') or die;

my $tmpdir = File::Temp::tempdir(CLEANUP => 1);
my $generator = Genome::Model::Command::Report::Generate->create(
    build => $build,
    report_name => 'Summary',
    directory => $tmpdir,
    force => 1,
);
ok($generator, 'create');
$generator->dump_status_messages(1);
ok($generator->execute, 'execute');

done_testing();
