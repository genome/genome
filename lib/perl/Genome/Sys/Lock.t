use strict;
use warnings;

use above 'Genome';
use Test::More tests => 4;
use List::Util qw(shuffle);

$ENV{GENOME_NESSY_SERVER} = 'http://nessy.gsc.wustl.edu/';

subtest 'optional backend lock does not cause failures' => sub {
    plan tests => 9;

    my $backends = [sort Genome::Sys::Lock->backends];
    is_deeply($backends, [qw(Genome::Sys::FileLock Genome::Sys::NessyLock)],
        'got expected backends');
    ok(Genome::Sys::Lock::is_mandatory('Genome::Sys::FileLock'), 'FileLock is mandatory');
    ok(!Genome::Sys::Lock::is_mandatory('Genome::Sys::NessyLock'), 'NessyLock is optional');

    no warnings 'redefine';
    local *Genome::Sys::NessyLock::lock = sub { 0 };
    use warnings 'redefine';

    my $resource_lock = 'Lock.t/' . random_string();
    my $lock = Genome::Sys::Lock->lock_resource(
        resource_lock => $resource_lock,
        block_sleep => 5,
        max_try => 0,
        wait_announce_interval => 10,
    );

    is($lock, $resource_lock, 'got lock even though NessyLock failed');
    ok(Genome::Sys::FileLock->has_lock($resource_lock), 'confirmed FileLock has the lock');
    ok(!Genome::Sys::NessyLock->has_lock($resource_lock), 'confirmed NessyLock does not have the lock');

    my $unlocked = Genome::Sys::Lock->unlock_resource(
        resource_lock => $lock,
    );
    ok($unlocked, 'unlocked');
    ok(!Genome::Sys::FileLock->has_lock($resource_lock), 'confirmed FileLock does not have the lock');
    ok(!Genome::Sys::NessyLock->has_lock($resource_lock), 'confirmed NessyLock does not have the lock');
};

subtest 'mandatory backend lock does cause failures' => sub {
    plan tests => 6;

    my $backends = [sort Genome::Sys::Lock->backends];
    is_deeply($backends, [qw(Genome::Sys::FileLock Genome::Sys::NessyLock)],
        'got expected backends');
    ok(Genome::Sys::Lock::is_mandatory('Genome::Sys::FileLock'), 'FileLock is mandatory');
    ok(!Genome::Sys::Lock::is_mandatory('Genome::Sys::NessyLock'), 'NessyLock is optional');

    no warnings 'redefine';
    local *Genome::Sys::FileLock::lock = sub { 0 };
    use warnings 'redefine';

    my $resource_lock = 'Lock.t/' . random_string();
    my $lock = Genome::Sys::Lock->lock_resource(
        resource_lock => $resource_lock,
        block_sleep => 5,
        max_try => 0,
        wait_announce_interval => 10,
    );

    is($lock, undef, 'failed to get lock');
    ok(!Genome::Sys::FileLock->has_lock($resource_lock), 'confirmed FileLock does not have the lock');
    ok(!Genome::Sys::NessyLock->has_lock($resource_lock), 'confirmed NessyLock does not have the lock');
};

subtest 'optional backend unlock does not cause failures' => sub {
    plan tests => 9;

    my $backends = [sort Genome::Sys::Lock->backends];
    is_deeply($backends, [qw(Genome::Sys::FileLock Genome::Sys::NessyLock)],
        'got expected backends');
    ok(Genome::Sys::Lock::is_mandatory('Genome::Sys::FileLock'), 'FileLock is mandatory');
    ok(!Genome::Sys::Lock::is_mandatory('Genome::Sys::NessyLock'), 'NessyLock is optional');

    my $resource_lock = 'Lock.t/' . random_string();
    my $lock = Genome::Sys::Lock->lock_resource(
        resource_lock => $resource_lock,
        block_sleep => 5,
        max_try => 0,
        wait_announce_interval => 10,
    );

    is($lock, $resource_lock, 'got lock');
    ok(Genome::Sys::FileLock->has_lock($resource_lock), 'confirmed FileLock has the lock');
    ok(Genome::Sys::NessyLock->has_lock($resource_lock), 'confirmed NessyLock has the lock');

    do {
        no warnings 'redefine';
        local *Genome::Sys::NessyLock::unlock = sub { 0 };
        use warnings 'redefine';

        my $unlocked = Genome::Sys::Lock->unlock_resource(
            resource_lock => $lock,
        );
        ok($unlocked, 'unlocked');
    };
    ok(!Genome::Sys::FileLock->has_lock($resource_lock), 'confirmed FileLock does not have the lock');
    ok(Genome::Sys::NessyLock->has_lock($resource_lock), 'confirmed NessyLock still has the lock');

    Genome::Sys::NessyLock->unlock($resource_lock);
};

subtest 'mandatory backend unlock does cause failures' => sub {
    plan tests => 9;

    my $backends = [sort Genome::Sys::Lock->backends];
    is_deeply($backends, [qw(Genome::Sys::FileLock Genome::Sys::NessyLock)],
        'got expected backends');
    ok(Genome::Sys::Lock::is_mandatory('Genome::Sys::FileLock'), 'FileLock is mandatory');
    ok(!Genome::Sys::Lock::is_mandatory('Genome::Sys::NessyLock'), 'NessyLock is optional');

    my $resource_lock = 'Lock.t/' . random_string();
    my $lock = Genome::Sys::Lock->lock_resource(
        resource_lock => $resource_lock,
        block_sleep => 5,
        max_try => 0,
        wait_announce_interval => 10,
    );

    is($lock, $resource_lock, 'got lock');
    ok(Genome::Sys::FileLock->has_lock($resource_lock), 'confirmed FileLock has the lock');
    ok(Genome::Sys::NessyLock->has_lock($resource_lock), 'confirmed NessyLock has the lock');

    do {
        no warnings 'redefine';
        local *Genome::Sys::FileLock::unlock = sub { 0 };
        use warnings 'redefine';

        my $unlocked = Genome::Sys::Lock->unlock_resource(
            resource_lock => $lock,
        );
        ok(!$unlocked, 'unlock failed');
    };
    ok(Genome::Sys::FileLock->has_lock($resource_lock), 'confirmed FileLock still has the lock');
    ok(!Genome::Sys::NessyLock->has_lock($resource_lock), 'confirmed NessyLock does not have the lock');

    Genome::Sys::FileLock->unlock(resource_lock => $resource_lock);
};

sub random_string {
    my @chars = map { (shuffle 'a'..'z', 'A'..'Z', 0..9)[0] } 1..10;
    return join('', @chars);
}
