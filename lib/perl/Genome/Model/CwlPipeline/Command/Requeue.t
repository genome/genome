#!/usr/bin/env genome-perl

BEGIN {
    $ENV{UR_DBI_NO_COMMIT} = 1;
};

use strict;
use warnings;

use above 'Genome';
use Test::More tests => 20;

use Genome::Test::Factory::Build;
use Genome::Test::Factory::Model::CwlPipeline;

my $class = 'Genome::Model::CwlPipeline::Command::Requeue';
use_ok($class) or die;

my $model = Genome::Test::Factory::Model::CwlPipeline->setup_object();

my $build = Genome::Test::Factory::Build->setup_object(model_id => $model->id, status => 'Running');
$build->run_by( Genome::Config::get('builder_run_as_user') );

my $user_build = Genome::Test::Factory::Build->setup_object(model_id => $model->id, status => 'Failed');
$user_build->run_by( 'not-' . Genome::Config::get('builder_run_as_user') );

run_test($build, undef, 'wrong status');
my @notes = $build->notes(header_text => 'Rebuild Requested');
is(scalar(@notes), 0, 'no note created');

$build->status('Failed');
run_test($build, undef, undef);
@notes = $build->notes(header_text => 'Rebuild Requested');
is(scalar(@notes), 1, 'created a note for the request');

my $reason = 'testing the requeue command';
run_test($build, $reason, undef);
@notes = $build->notes(header_text => 'Rebuild Requested');
is(scalar(@notes), 2, 'created a second note given the second request');
my @with_reason = grep { $_->body_text eq $reason } @notes;
is(scalar(@with_reason), 1, 'found note with supplied reason');

ok(UR::Context->commit(), 'can sync rebuild request to db');
$build->rebuild_requested(0);
ok(UR::Context->commit(), 'can replace value later (such as when restarting build');

run_test($user_build, $reason, 'wrong user');
my @user_notes = $user_build->notes(header_text => 'Rebuild Requested');
is(scalar(@user_notes), 0, 'no note created');

done_testing();

sub run_test {
    my ($build, $reason, $failure_reason) = @_;

    my @params = ( builds => [$build] );
    if ($reason) {
        push @params, reason => $reason;
    }

    my $cmd = $class->create(@params);
    isa_ok($cmd, $class, 'created command');
    ok($cmd->execute, 'executed command');

    if ($failure_reason) {
        ok(!$build->rebuild_requested, 'rebuild was not requested -- ' . $failure_reason);
    } else {
        ok($build->rebuild_requested, 'rebuild was requested');
    }

    Genome::Sys::Lock->release_all(); #setting metrics creates a lock; normally the command would commit but not here in testing
}

