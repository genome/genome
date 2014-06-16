#!/usr/bin/env genome-perl

use strict;
use warnings;

use above 'Genome';

use Data::Dumper;
use Genome::Model::Build::MetagenomicComposition16s::TestBuildFactory;
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
my ($build, $example_build) = Genome::Model::Build::MetagenomicComposition16s::TestBuildFactory->build_with_example_build_for_454;
my $model = $build->model;

#< FAIL >#
# fail - no models or builds
my $cmd = Genome::Model::MetagenomicComposition16s::Command::Tester->create();
ok($cmd, 'create');
$cmd->dump_status_messages(1);
ok(!$cmd->execute, 'execute failed w/o builds');

is($build->status('Succeeded'), 'Succeeded', 'build is succeeded');

# execute ok - models and builds
$cmd = Genome::Model::MetagenomicComposition16s::Command::Tester->create(
    models => [$model],
);
ok($cmd, 'create');
$cmd->dump_status_messages(1);
ok($cmd->execute, 'execute');
is_deeply([$cmd->_builds], [$build], 'builds from cmd') 
    or diag(
        Data::Dumper::Dumper([$cmd->_builds],[$build])
    );

$cmd = Genome::Model::MetagenomicComposition16s::Command::Tester->create(
    builds => [$build],
);
ok($cmd, 'create');
$cmd->dump_status_messages(1);
ok($cmd->execute, 'execute');
is_deeply([$cmd->_builds], [$build], 'builds from cmd');

done_testing();
