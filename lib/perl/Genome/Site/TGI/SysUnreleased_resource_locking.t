#!/usr/bin/env genome-perl
use strict;
use warnings;
use Test::More tests => 5;
use above 'Genome';

my $tmp_dir = Genome::Sys->create_temp_directory();
ok($tmp_dir, "created temp dir ($tmp_dir)");

my @common_params = (lock_directory => $tmp_dir, resource_id => "foo", block_sleep => 0);

$SIG{CHLD} = sub { wait };
my $child_pid = fork;
if ($child_pid == 0) { # child thread
    print "CHILD: Locking $tmp_dir/foo...\n";
    my $child_lock = Genome::Sys->lock_resource(@common_params, max_try => 0);
    print "CHILD: Sleeping for two seconds...\n";
    sleep(2);
    print "CHILD: Unlocking $tmp_dir/foo...\n";
    Genome::Sys->unlock_resource(resource_lock => $child_lock);
    print "CHILD: Exiting...\n";
    exit 0;
}
else { # parent thread
    sleep(1);
    print "PARENT: Trying to lock $tmp_dir/foo...\n";
    my $parent_lock = Genome::Sys->lock_resource(@common_params, max_try => 0);
    is($parent_lock, undef, 'correctly failed to get lock on temp dir while child process has it locked');

    waitpid($child_pid, 0);
}

ok(Genome::Sys->lock_resource(@common_params, max_try => 2), 'locked temp dir once child process finished');
ok(Genome::Sys->lock_resource(@common_params, max_try => 2), 'locked temp dir even though I already locked it');
my ($last_warning) = Genome::Sys->warning_message;
is($last_warning, "Looks like I'm waiting on my own lock, forcing unlock...", 'got warning about waiting on own lock');

