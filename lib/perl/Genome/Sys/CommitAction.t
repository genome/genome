#!/usr/bin/env genome-perl

use strict;
use warnings;

use above 'Genome';
use Test::More tests => 2;
use Test::Exception;

subtest basic => sub {
    plan tests => 7;

    my($sync_called, $commit_called) = (0, 0);

    my $action = Genome::Sys::CommitAction->create(
                    on_sync => sub { $sync_called++ },
                    on_commit => sub { $commit_called++ },
                );
    ok($action, 'Created Genome::Sys::CommitAction');

    ok(UR::Context->commit(), 'commit');

    is($sync_called, 1, 'on_sync callback run');
    is($commit_called, 1, 'on_commit_callback run');

    ok(UR::Context->commit(), 'commit again');

    is($sync_called, 1, 'on_sync callback was not run again');
    is($commit_called, 1, 'on_commit_callback was not run again');
};

subtest rollback => sub {
    plan tests => 5;

    my($sync_called, $commit_called, $sync_failed_called) = (0,0,0);

    my $action = Genome::Sys::CommitAction->create(
                    on_sync => sub { $sync_called++; die "in on_sync\n"; },
                    on_commit => sub { $commit_called++ },
                    on_sync_fail => sub { $sync_failed_called++ },
                );
    ok($action, 'Created Genome::Sys::CommitAction');

    throws_ok { UR::Context->commit() }
              qr/in on_sync/,
              'exception in on_sync';

    is($sync_called, 1, 'on_sync callback not run');
    is($commit_called, 0, 'on_commit callback not run');
    is($sync_failed_called, 1, 'on_sync_fail callback was run');
};
