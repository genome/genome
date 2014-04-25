use strict;
use warnings;

BEGIN {
    $ENV{GENOME_NESSY_SERVER} = 'http://nessy.gsc.wustl.edu/';
};

use above 'Genome';
use Test::More tests => 8;
use List::Util qw(shuffle);

my @backends = sort Genome::Sys::Lock->backends;
is(scalar(@backends), 2, 'got two backends');

my ($filelock) = grep { $_ eq 'Genome::Sys::FileLock' } @backends;
ok($filelock->is_mandatory(), 'FileLock is mandatory');

my ($nessylock) = grep { $_->can('blessed') } @backends;
ok(!$nessylock->is_mandatory(), 'NessyBackend is optional');

subtest 'optional backend lock does not cause failures' => sub {
    plan tests => 6;

    no warnings 'redefine';
    local *Genome::Sys::Lock::NessyBackend::lock = sub { 0 };
    use warnings 'redefine';

    my $resource_lock = 'Lock.t/' . random_string();
    my $lock = Genome::Sys::Lock->lock_resource(
        resource_lock => $resource_lock,
        block_sleep => 5,
        max_try => 0,
        wait_announce_interval => 10,
    );

    is($lock, $resource_lock, 'got lock even though NessyBackend failed');
    ok(Genome::Sys::FileLock->has_lock($resource_lock), 'confirmed FileLock has the lock');
    ok(!$nessylock->has_lock($resource_lock), 'confirmed NessyBackend does not have the lock');

    my $unlocked = Genome::Sys::Lock->unlock_resource(
        resource_lock => $lock,
    );
    ok($unlocked, 'unlocked');
    ok(!Genome::Sys::FileLock->has_lock($resource_lock), 'confirmed FileLock does not have the lock');
    ok(!$nessylock->has_lock($resource_lock), 'confirmed NessyBackend does not have the lock');
};

subtest 'mandatory backend lock does cause failures' => sub {
    plan tests => 3;

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
    ok(!$nessylock->has_lock($resource_lock), 'confirmed NessyBackend does not have the lock');
};

subtest 'optional backend unlock does not cause failures' => sub {
    plan tests => 6;

    my $resource_lock = 'Lock.t/' . random_string();
    my $lock = Genome::Sys::Lock->lock_resource(
        resource_lock => $resource_lock,
        block_sleep => 5,
        max_try => 0,
        wait_announce_interval => 10,
    );

    is($lock, $resource_lock, 'got lock');
    ok(Genome::Sys::FileLock->has_lock($resource_lock), 'confirmed FileLock has the lock');
    ok($nessylock->has_lock($resource_lock), 'confirmed NessyBackend has the lock');

    do {
        no warnings 'redefine';
        local *Genome::Sys::Lock::NessyBackend::unlock = sub { 0 };
        use warnings 'redefine';

        my $unlocked = Genome::Sys::Lock->unlock_resource(
            resource_lock => $lock,
        );
        ok($unlocked, 'unlocked');
    };
    ok(!Genome::Sys::FileLock->has_lock($resource_lock), 'confirmed FileLock does not have the lock');
    ok($nessylock->has_lock($resource_lock), 'confirmed NessyBackend still has the lock');

    $nessylock->unlock($resource_lock);
};

subtest 'mandatory backend unlock does cause failures' => sub {
    plan tests => 6;

    my $resource_lock = 'Lock.t/' . random_string();
    my $lock = Genome::Sys::Lock->lock_resource(
        resource_lock => $resource_lock,
        block_sleep => 5,
        max_try => 0,
        wait_announce_interval => 10,
    );

    is($lock, $resource_lock, 'got lock');
    ok(Genome::Sys::FileLock->has_lock($resource_lock), 'confirmed FileLock has the lock');
    ok($nessylock->has_lock($resource_lock), 'confirmed NessyBackend has the lock');

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
    ok(!$nessylock->has_lock($resource_lock), 'confirmed NessyBackend does not have the lock');

    Genome::Sys::FileLock->unlock(resource_lock => $resource_lock);
};

subtest 'mandatory twin NessyBackends' => sub {
    plan tests => 2;

    my @backends = (
        Genome::Sys::Lock::NessyBackend->new(
            url => 'http://nessy.gsc.wustl.edu/',
            is_mandatory => 1,
        ),
        Genome::Sys::Lock::NessyBackend->new(
            url => 'http://nessy.gsc.wustl.edu/',
            is_mandatory => 1,
        ),
    );
    no warnings 'redefine';
    local *Genome::Sys::Lock::backends = sub { @backends };
    use warnings 'redefine';

    my $resource_lock = 'Lock.t/' . random_string();
    my $lock = Genome::Sys::Lock->lock_resource(
        resource_lock => $resource_lock,
        block_sleep => 5,
        max_try => 0,
        wait_announce_interval => 10,
    );
    is($lock, undef, 'failed to get lock');
    ok(!(grep { $_->has_lock($resource_lock) } @backends), 'neither backend locked the resource');
};

sub random_string {
    my @chars = map { (shuffle 'a'..'z', 'A'..'Z', 0..9)[0] } 1..10;
    return join('', @chars);
}
