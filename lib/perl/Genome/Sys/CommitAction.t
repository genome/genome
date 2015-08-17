#!/usr/bin/env genome-perl

use strict;
use warnings;

use above 'Genome';
use Test::More tests => 2;
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

subtest rollback => sub {
    plan tests => 8;

    my($sync_called, $commit_called, $sync_failed_called) = (0,0,0);

    my $action1 = Genome::Sys::CommitAction->create(
                    on_sync => sub { $sync_called++; die "in on_sync\n"; },
                    on_commit => sub { $commit_called++ },
                    on_sync_fail => sub { $sync_failed_called++ },
                );
    ok($action1, 'Created Genome::Sys::CommitAction');

    my $action2 = Genome::Sys::CommitAction->create(
                    on_sync => sub { 1 },
                  );
    ok($action2, 'Created second CommitAction');

    throws_ok { UR::Context->commit() }
              qr/in on_sync/,
              'exception in on_sync';

    is($sync_called, 1, 'on_sync callback not run');
    is($commit_called, 0, 'on_commit callback not run');
    is($sync_failed_called, 1, 'on_sync_fail callback was run');

    isa_ok($action1, 'Genome::Sys::CommitAction', 'First action is a Genome::Sys::CommitAction');
    isa_ok($action2, 'Genome::Sys::CommitAction', 'Second action is a Genome::Sys::CommitAction');
};
