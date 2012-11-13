#!/usr/bin/env genome-perl

use strict;
use warnings;

use above 'Genome';

use Data::Dumper;
use Genome::Model::MetagenomicComposition16s::Test;
use Test::More;

use_ok('Genome::Model::MetagenomicComposition16s::Command');

# fake class and execute to test base class
class Genome::Model::MetagenomicComposition16s::Command::Tester {
    is => 'Genome::Model::MetagenomicComposition16s::Command',
};
sub Genome::Model::MetagenomicComposition16s::Command::Tester::execute { 
    my $self = shift;
    return $self->_builds;
}

# model
my $model = Genome::Model::MetagenomicComposition16s::Test->model_for_sanger;
ok($model, 'Got mock MC16s sanger model');
my $build = Genome::Model::Build::MetagenomicComposition16s->create(
    model => $model,
);
ok($build, 'Added build to model');
ok($build->the_master_event->date_completed(UR::Context->current->now), 'build has date completed');

#< FAIL >#
# fail - no builds
my $cmd = Genome::Model::MetagenomicComposition16s::Command::Tester->create();
ok($cmd, 'create');
$cmd->dump_status_messages(1);
ok(!$cmd->execute, 'execute failed w/o builds');

# fail - model doesn't have a build
$cmd = Genome::Model::MetagenomicComposition16s::Command::Tester->create(
    models => [$model],
);
ok($cmd, 'create');
$cmd->dump_status_messages(1);
ok(!$cmd->execute, 'execute failed w/ model w/o builds');

is($build->the_master_event->event_status('Succeeded'), 'Succeeded', 'build is succeeded');

# execute ok - models and builds
$cmd = Genome::Model::MetagenomicComposition16s::Command::Tester->create(
    models => [$model],
);
ok($cmd, 'create');
$cmd->dump_status_messages(1);
ok($cmd->execute, 'execute');
is_deeply([$cmd->_builds], [$build], 'builds from cmd');

$cmd = Genome::Model::MetagenomicComposition16s::Command::Tester->create(
    builds => [$build],
);
ok($cmd, 'create');
$cmd->dump_status_messages(1);
ok($cmd->execute, 'execute');
is_deeply([$cmd->_builds], [$build], 'builds from cmd');

done_testing();
exit;

