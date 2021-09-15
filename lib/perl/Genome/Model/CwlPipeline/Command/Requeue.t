#!/usr/bin/env genome-perl

BEGIN {
    $ENV{UR_DBI_NO_COMMIT} = 1;
};

use strict;
use warnings;

use Sub::Override;

use above 'Genome';
use Test::More tests => 36;

use Genome::Test::Factory::AnalysisProject;
use Genome::Test::Factory::Build;
use Genome::Test::Factory::Model::CwlPipeline;

my $class = 'Genome::Model::CwlPipeline::Command::Requeue';
use_ok($class) or die;

my $anp = Genome::Test::Factory::AnalysisProject->setup_object();
$anp->created_by('-nonexistent-user-');
my $model = Genome::Test::Factory::Model::CwlPipeline->setup_object();

Genome::Config::AnalysisProject::ModelBridge->create(analysis_project_id => $anp->id, model_id => $model->id);

my $build = Genome::Test::Factory::Build->setup_object(model_id => $model->id, status => 'Running');
$build->run_by( Genome::Config::get('builder_run_as_user') );

my $user_build = Genome::Test::Factory::Build->setup_object(model_id => $model->id, status => 'Failed');
$user_build->run_by( 'not-' . Genome::Config::get('builder_run_as_user') );

run_test($build, undef, 'wrong status');
my @notes = $build->notes(header_text => 'Action Requested: restart');
is(scalar(@notes), 0, 'no note created');

$build->status('Failed');
run_test($build, undef, undef);
@notes = $build->notes(header_text => 'Action Requested: restart');
is(scalar(@notes), 1, 'created a note for the request');

my $reason = 'testing the requeue command';
run_test($build, $reason, undef);
@notes = $build->notes(header_text => 'Action Requested: restart');
is(scalar(@notes), 2, 'created a second note given the second request');
my @with_reason = grep { $_->body_text eq $reason } @notes;
is(scalar(@with_reason), 1, 'found note with supplied reason');

ok(UR::Context->commit(), 'can sync rebuild request to db');
$build->action_requested(0);
ok(UR::Context->commit(), 'can replace value later (such as when restarting build');

run_test($user_build, $reason, 'wrong user');
my @user_notes = $user_build->notes(header_text => 'Action Requested: restart');
is(scalar(@user_notes), 0, 'no note created');

run_test($user_build, $reason, 'wrong user', 'stop');
@user_notes = $user_build->notes(header_text => 'Action Requested: stop');
is(scalar(@user_notes), 0, 'no note created');

use Genome::Sys;
my $guard = Sub::Override->new('Genome::Sys::current_user_is_admin', sub { return 0; });

run_test($build, undef, 'wrong AnP owner', 'stop');
@notes = $build->notes(header_text => 'Action Requested: stop');
is(scalar(@notes), 0, 'no note created');

$anp->created_by(Genome::Sys->username);

run_test($build, undef, undef, 'stop');
@notes = $build->notes(header_text => 'Action Requested: stop');
is(scalar(@notes), 1, 'created a note for the request');

run_test($build, undef, undef, '0');
@notes = $build->notes(header_text => 'Action Requested: 0');
is(scalar(@notes), 1, 'created a note for cancelling the request');

undef $guard;

done_testing();

sub run_test {
    my ($build, $reason, $failure_reason, $action) = @_;

    my @params = ( builds => [$build] );
    if (defined $reason) {
        push @params, reason => $reason;
    }
    if (defined $action) {
        push @params, action => $action;
    }

    my $cmd = $class->create(@params);
    isa_ok($cmd, $class, 'created command');
    ok($cmd->execute, 'executed command');

    if ($failure_reason) {
        ok(!$build->action_requested, 'rebuild was not requested -- ' . $failure_reason);
    } else {
        is($build->action_requested, ($action // 'restart'), 'action was requested');
    }

    Genome::Sys::Lock->release_all(); #setting metrics creates a lock; normally the command would commit but not here in testing
}

