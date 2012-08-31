#!/usr/bin/env genome-perl

use strict;
use warnings;

use above 'Genome';

use Test::More;

my $build = Genome::Model::Build->get(107664200); # build for apipe-test-03-MC16s
ok($build, 'Got MC16s build') or die;

my $generator = Genome::Model::Report::BuildInitialized->create(build_id => $build->id);
ok($generator, 'create');
my $report = $generator->generate_report;
ok($report, 'generate report');

done_testing();
exit;

