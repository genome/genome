#!/usr/bin/env genome-perl

# Testing a theory about LSF state latency.
sleep 5;

use strict;
use warnings;

$Genome::Sys::IS_TESTING=1;
BEGIN {
    $ENV{UR_DBI_NO_COMMIT} = 1;
};

use above 'Genome';
use Carp qw(croak);
use Data::Dumper;
use File::Path;
use File::Slurp;
use File::Spec;
use File::Temp;
use IO::Handle qw();
use MIME::Base64;
use POSIX ":sys_wait_h";
use Socket qw(AF_UNIX SOCK_STREAM PF_UNSPEC);
use Test::More;
use Time::HiRes qw(gettimeofday);

require_ok('Genome::Sys::Lock');

my $tmp_dir = File::Temp::tempdir('Genome-Utility-FileSystem-writetest-XXXXX', CLEANUP => 1, TMPDIR => 1);

my $bogus_id = '-55555';
$tmp_dir = File::Temp::tempdir(CLEANUP => 1);
my $sub_dir = $tmp_dir .'/sub/dir/test';
ok(! -e $sub_dir,$sub_dir .' does not exist');
ok(Genome::Sys->create_directory($sub_dir),'create directory');
ok(-d $sub_dir,$sub_dir .' is a directory');

test_locking(successful => 1,
    message => 'lock resource_id '. $bogus_id,
    lock_directory => $tmp_dir,
    resource_id => $bogus_id,);

{
    # Override is_my_lock_target so this doesn't carp for trying to lock same
    # resource twice.  This replaces the old option to wait_on_self.
    no warnings 'redefine';
    local *Genome::Sys::Lock::is_my_lock_target = sub { 0 };
    test_locking(successful => 0,
        message => 'failed lock resource_id '. $bogus_id,
        lock_directory => $tmp_dir,
        resource_id => $bogus_id,
        max_try => 1,
        block_sleep => 3,);
}

ok(Genome::Sys->unlock_resource(
        lock_directory => $tmp_dir,
        resource_id => $bogus_id,
    ), 'unlock resource_id '. $bogus_id);
my $init_lsf_job_id = $ENV{'LSB_JOBID'};
{
    local $ENV{'LSB_JOBID'};
    $ENV{'LSB_JOBID'} = 1;
    test_locking(successful => 1,
        message => 'lock resource with bogus lsf_job_id',
        lock_directory => $tmp_dir,
        resource_id => $bogus_id,);
    {
        # Override is_my_lock_target so this doesn't carp for trying to lock same
        # resource twice.  This replaces the old option to wait_on_self.
        no warnings 'redefine';
        local *Genome::Sys::Lock::is_my_lock_target = sub { 0 };
        test_locking(
            successful=> 1,
            message => 'lock resource with removing invalid lock with bogus lsf_job_id first',
            lock_directory => $tmp_dir,
            resource_id => $bogus_id,
            max_try => 1,
            block_sleep => 3,);
    }
    ok(Genome::Sys->unlock_resource(
            lock_directory => $tmp_dir,
            resource_id => $bogus_id,
        ), 'unlock resource_id '. $bogus_id);
# TODO: add skip test but if we are on a blade, lets see that the locking works correctly
# Above the test is that old bogus locks can get removed when the lsf_job_id no longer exists
# We should test that while an lsf_job_id does exist (ie. our current job) we still hold the lock
    SKIP: {
        skip 'only test the state of the lsf job if we are running on a blade with a job id',
        3 unless ($init_lsf_job_id);
        $ENV{'LSB_JOBID'} = $init_lsf_job_id;
        ok(Genome::Sys->lock_resource(
                lock_directory => $tmp_dir,
                resource_id => $bogus_id,
            ),'lock resource with real lsf_job_id');
        {
            # Override is_my_lock_target so this doesn't carp for trying to lock same
            # resource twice.  This replaces the old option to wait_on_self.
            no warnings 'redefine';
            local *Genome::Sys::Lock::is_my_lock_target = sub { 0 };
            ok(!Genome::Sys->lock_resource(
                    lock_directory => $tmp_dir,
                    resource_id => $bogus_id,
                    max_try => 1,
                    block_sleep => 3,
                ),'failed lock resource with real lsf_job_id blocking');
        }
        ok(Genome::Sys->unlock_resource(
                lock_directory => $tmp_dir,
                resource_id => $bogus_id,
            ), 'unlock resource_id '. $bogus_id);
    };
}

# RACE CONDITION
my $base_dir = File::Temp::tempdir("Genome-Utility-FileSystem-RaceCondition-XXXX", CLEANUP=>1, TMPDIR => 1);

my $resource = "/tmp/Genome-Utility-Filesystem.test.resource.$$";

if (defined $ENV{'LSB_JOBID'} && $ENV{'LSB_JOBID'} eq "1") {
    delete $ENV{'LSB_JOBID'}; 
}

my @pids;

my $children = 20;

for my $child (1...$children) {
    my $pid;
    if ($pid = UR::Context::Process->fork()) {
        push @pids, $pid;
    } else {
        my $output_offset = $child;
        my $tempdir = $base_dir;
        my $output_file = $tempdir . "/" . $output_offset;

        my $fh = new IO::File(">>$output_file");

        my $lock = Genome::Sys->lock_resource(
            resource_lock => $resource,
            block_sleep   => 1,
            wait_announce_interval => 30,
        );
        unless ($lock) {
            print_event($fh, "LOCK_FAIL", "Failed to get a lock" );
            $fh->close;
            exit(1);
        }

        #sleep for half a second before printing this, to let a prior process catch up with
        #printing its "unlocked" message.  sometimes we get a lock (properly) in between the time
        #someone else has given up the lock, but before it had a chance to report that it did.
        select(undef, undef, undef, 0.50);
        print_event($fh, "LOCK_SUCCESS", "Successfully got a lock" );
        sleep 2;

        unless (Genome::Sys->unlock_resource(resource_lock => $resource)) {
            print_event($fh, "UNLOCK_FAIL", "Failed to release a lock" );
            $fh->close;
            exit(1);
        }
        print_event($fh, "UNLOCK_SUCCESS", "Successfully released my lock" );

        $fh->close;
        exit(0);
    }
}

for my $pid (@pids) {
    my $status = waitpid $pid, 0;
}

my @event_log;

for my $child (1...$children) {
    my $report_log = "$base_dir/$child";
    ok (-e $report_log, "Expected to see a report for $report_log");
    if (!-e $report_log) {
        die "Expected output file did not exist";
    } else {
        my @lines = read_file($report_log);
        for (@lines) {
            my ($time, $event, $pid, $msg) = split /\t/;
            my $report = {etime => $time, event=>$event, pid=>$pid, msg=>$msg};
            push @event_log, $report;
        } 
    }
}

ok(scalar @event_log  == 2*$children, "read in got 2 lock/unlock events for each child.");

@event_log = sort {$a->{etime} <=> $b->{etime}} @event_log;
my $last_event;
for (@event_log) {
    if (defined $last_event) {
        my $valid_next_event = ($last_event eq "UNLOCK_SUCCESS" ? "LOCK_SUCCESS" : "UNLOCK_SUCCESS");
        ok($_->{event} eq $valid_next_event, sprintf("Last lock event was a %s so we expected a to see %s, got a %s", $last_event, $valid_next_event, $_->{event}));
    }
    $last_event = $_->{event};
    printf("%s\t%s\t%s\n", $_->{etime}, $_->{pid}, $_->{event});
}

my $tmp_dir2 = Genome::Sys->create_temp_directory();
ok($tmp_dir2, "created temp dir ($tmp_dir2)");

my @common_params = (lock_directory => $tmp_dir2, resource_id => "foo", block_sleep => 0);

my $child_sleep = Genome::Sys::Lock->min_timeout() + 2;
my $child_pid = UR::Context::Process->fork;
if ($child_pid == 0) { # child thread
    print "CHILD: Locking $tmp_dir2/foo...\n";
    my $child_lock = Genome::Sys->lock_resource(@common_params, max_try => 0);
    unless ($child_lock) {
        Carp::croak('CHILD: Failed to get lock.');
    }
    print "CHILD: Sleeping for two seconds...\n";
    sleep($child_sleep);
    print "CHILD: Unlocking $tmp_dir2/foo...\n";
    Genome::Sys->unlock_resource(resource_lock => $child_lock);
    print "CHILD: Exiting...\n";
    exit 0;
}
else { # parent thread
    sleep(1);
    print "PARENT: Trying to lock $tmp_dir2/foo...\n";
    my $parent_lock = Genome::Sys->lock_resource(@common_params, max_try => 0);
    is($parent_lock, undef, 'correctly failed to get lock on temp dir while child process has it locked');

    waitpid($child_pid, 0);
}

ok(Genome::Sys->lock_resource(@common_params, max_try => 2), 'locked temp dir once child process finished');

done_testing();

###

sub test_locking {
    my %params = @_;
    my $successful = delete $params{successful};
    die unless defined($successful);
    my $message = delete $params{message};
    die unless defined($message);

    my $lock = Genome::Sys->lock_resource(%params);
    if ($successful) {
        ok($lock,$message);
        if ($lock) { return $lock; }
    } else {
        ok(!$lock,$message);
    }
    return;
}

sub print_event {
    my $fh = shift;
    my $info = shift;
    my $msg  = shift;

    my ( $seconds, $ms ) = gettimeofday();
    $ms = sprintf("%06d",$ms);
    my $time = "$seconds.$ms";

    my $tp = sprintf( "%s\t%s\t%s\t%s", $time, $info, $$, $msg );

    print $fh $tp, "\n";
    print $tp, "\n";
}

sub fork_lock_resource {
    my %args = @_;

    socketpair(my $CHILD, my $PARENT, AF_UNIX, SOCK_STREAM, PF_UNSPEC)
        or die "socketpair: $!";

    $CHILD->autoflush(1);
    $PARENT->autoflush(1);

    my $pid = UR::Context::Process->fork();
    if ($pid) {
        close $PARENT;
        waitfor($CHILD, 'READY');
        $CHILD->say('LOCK');
        my $lock = getline($CHILD);
        $CHILD->say('ACK');

        close $CHILD;
        waitpid($pid, 0);

        return $lock;
    } else {
        close $CHILD;
        Genome::Sys::Lock->clear_state();
        $PARENT->say('READY');

        waitfor($PARENT, 'LOCK');
        my $lock = Genome::Sys::Lock->lock_resource(%args);
        $PARENT->say($lock);

        waitfor($PARENT, 'ACK');
        if ($lock) {
            Genome::Sys::Lock->unlock_resource(resource_lock => $lock);
        }
        close $PARENT;
        exit(0);
    }
}

sub getline {
    my $handle = shift;
    unless ($handle) {
        croak('handle required');
    }
    my $line = <$handle>;
    chomp $line;
    return $line;
}

sub waitfor {
    my ($handle, $message) = @_;
    unless ($handle) {
        croak('handle required');
    }
    unless (defined $message) {
        croak('message required');
    }
    my $line = getline($handle);
    unless ($line) {
        croak("empty message");
    }
    unless ($line eq $message) {
        croak("unexpected message: $line");
    }
    return 1;
}
