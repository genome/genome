#!/usr/bin/env genome-perl

use strict;
use warnings;

use above 'Genome';

use Test::More;

use_ok('Genome::Model::MetagenomicComposition16s::Command::ListRuns') or die;

# model/build
my $model = Genome::Model::MetagenomicComposition16s->get(name => 'apipe-test-mc16s-sanger');
ok($model, 'Got mock mc16s sanger model');

my $tmpdir = File::Temp::tempdir(CLEANUP => 1);

# ok - list w/ model name
my $cmd = Genome::Model::MetagenomicComposition16s::Command::ListRuns->create(
    models => [$model],
);
ok($cmd, 'create list runs');
$cmd->dump_status_messages(1);
ok($cmd->execute, 'Execute list runs ok');

done_testing();
