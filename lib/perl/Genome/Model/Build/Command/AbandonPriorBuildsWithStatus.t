#!/usr/bin/env genome-perl

use strict;
use warnings;

BEGIN {
    $ENV{UR_DBI_NO_COMMIT} = 1;
};

use above 'Genome';

use Genome::Test::Factory::Model::CwlPipeline;
use Genome::Test::Factory::Build;

use Test::More tests => 7;

my $pkg = 'Genome::Model::Build::Command::AbandonPriorBuildsWithStatus';

use_ok($pkg);

my $model = Genome::Test::Factory::Model::CwlPipeline->setup_object;
for (('Failed')x3, 'Succeeded', 'Failed', ('Succeeded')x2, 'Scheduled', ('Succeeded')x3) {
    my $build = Genome::Test::Factory::Build->setup_object(
        model_id => $model->id,
        status => $_,
    );
    sleep 2; #our build timestamps only have 1s resolution so let's not make them all at the same time
}
is(scalar(@{[$model->builds]}), 11, 'created 11 test builds');

my $last_complete_build = $model->last_complete_build;

my $cmd = $pkg->create(
    models => [$model],
    status => 'Succeeded'
);
isa_ok($cmd, $pkg, 'created command');
ok($cmd->execute, 'executed command');

my %statuses;
for my $build ($model->builds) {
    $statuses{$build->status}++;
}

is($statuses{Succeeded}, 1, 'left one Succeeded build');
is($statuses{Abandoned}, 5, 'abandoned other Succeeded builds');
is($last_complete_build->status, 'Succeeded', 'the remaining Succeeded build is the most recent one');
