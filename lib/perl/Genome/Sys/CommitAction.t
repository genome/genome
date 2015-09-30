#!/usr/bin/env genome-perl

use strict;
use warnings;

use above 'Genome';
use Test::More tests => 3;
use Test::Exception;

subtest basic => sub {
    plan tests => 11;

    my($sync_called, $commit_called, $sync_received_data, $commit_received_data) = (0, 0, undef);
    my $passed_data = 'hi there';

    my $action = Genome::Sys::CommitAction->create(
                    on_sync => sub { $sync_called++; $sync_received_data = shift; },
                    on_commit => sub { $commit_called++; $commit_received_data = shift; },
                    data => $passed_data,
                );
    ok($action, 'Created Genome::Sys::CommitAction');

    ok(UR::Context->commit(), 'commit');

    is($sync_called, 1, 'on_sync callback run');
    is($commit_called, 1, 'on_commit_callback run');
    is($sync_received_data, $passed_data, 'callback got data');
    is($commit_received_data, $passed_data, 'callback got data');

    ok(UR::Context->commit(), 'commit again');

    is($sync_called, 1, 'on_sync callback was not run again');
    is($commit_called, 1, 'on_commit_callback was not run again');

    isa_ok($action, 'UR::DeletedRef', 'CommitAction is deleted');
    my @actions = Genome::Sys::CommitAction->is_loaded();
    is(scalar(@actions), 0, 'No CommitActions left in memory');
};

subtest 'lifecycle' => sub {
    plan tests => 15;

    my($sync_called, $commit_called, $sync_failed_called, $rollback_called) = (0,0,0,0);

    my $action1 = Genome::Sys::CommitAction->create(
                    on_sync => sub { $sync_called++; die "in on_sync\n"; },
                    on_commit => sub { $commit_called++ },
                    on_sync_fail => sub { $sync_failed_called++ },
                    on_rollback => sub { $rollback_called++ },
                );
    ok($action1, 'Created Genome::Sys::CommitAction');

    my $action2 = Genome::Sys::CommitAction->create(
                    on_sync => sub { $sync_called++ },
                    on_commit => sub { $commit_called++ },
                    on_sync_fail => sub { $sync_failed_called++ },
                    on_rollback => sub { $rollback_called++ },
                  );
    ok($action2, 'Created second CommitAction');

    throws_ok { UR::Context->commit() }
              qr/in on_sync/,
              'exception in on_sync';

    is($commit_called, 0, 'on_commit callback not run');
    ok($sync_failed_called > 1, 'on_sync_fail callback was run');
    is($rollback_called, 0, 'on_rollback was not run');

    isa_ok($action1, 'Genome::Sys::CommitAction', 'First action');
    isa_ok($action2, 'Genome::Sys::CommitAction', 'Second action');


    ($sync_called, $commit_called, $sync_failed_called, $rollback_called) = (0,0,0,0);
    $action1->on_sync( sub { $sync_called++ });
    ok(UR::Context->commit(), 'commit');
    is($sync_called, 2, 'on_sync called twice');
    is($commit_called, 2, 'on_commit called twice');
    is($sync_failed_called, 0, 'on_sync_failed was not run');
    is($rollback_called, 0, 'on_rollback was not run');
    isa_ok($action1, 'UR::DeletedRef', 'First action');
    isa_ok($action2, 'UR::DeletedRef', 'Second action');
};

subtest 'rollback' => sub {
    plan tests => 4;

    my($sync_fail, $rollback) = (0,0);
    my $action = Genome::Sys::CommitAction->create(
                    on_sync_fail => sub { $sync_fail++ },
                    on_rollback  => sub { $rollback++ },
                );

    ok(UR::Context->rollback, 'rollback');
    is($rollback, 1, 'on_rollback called');
    is($sync_fail, 0, 'on_sync_fail not called');
    isa_ok($action, 'UR::DeletedRef');
};
