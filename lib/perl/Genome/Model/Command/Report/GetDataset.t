#!/usr/bin/env genome-perl

use strict;
use warnings;

use above 'Genome';

use Test::More;

use_ok('Genome::Model::Command::Report::GetDataset') or die;

my $build = Genome::Model::Build->get(107664200); # build for apipe-test-03-MC16s
ok($build, 'Got MC16s build') or die;

print $build->reports_directory."\n";
my $get_ds = Genome::Model::Command::Report::GetDataset->create(
    build => $build,
    report_name => 'Summary',
    dataset_name => 'stats',
    output_type => 'csv',
);
ok($get_ds, 'create');
$get_ds->dump_status_messages(1);
ok($get_ds->execute, 'execute');

done_testing();
