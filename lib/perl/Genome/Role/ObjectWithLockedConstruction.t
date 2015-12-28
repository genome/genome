use strict;
use warnings;

use Genome;

use Test::More tests => 3;
use Test::Exception;

subtest 'required methods' => sub {
    plan tests => 2;

    throws_ok { class RequiredMethodsFail {
                    roles => 'Genome::Role::ObjectWithLockedConstruction',
                } }
        qr/missing required property or method 'lock_id'/,
        'Role requires lock_id method';

    sub RequiredMethodsSucceed::lock_id { 1 }
    lives_ok { class RequiredMethodsSucceed {
                    roles => 'Genome::Role::ObjectWithLockedConstruction',
            } }
            'Compose role';
};

subtest 'create with no commit on' => sub {
    plan tests => 4;

    local $ENV{UR_DBI_NO_COMMIT} = 1;
    my $lock_id_called = 0;

    my $locker = LockHijacker->new();

    { no warnings 'once';
        *TestCreateNoCommit::lock_id = sub { ++$lock_id_called };
    }
    class TestCreateNoCommit {
        roles => ['Genome::Role::ObjectWithLockedConstruction'],
    };

    my $o = TestCreateNoCommit->create();
    ok($o, 'Created object');
    is($lock_id_called, 0, 'method lock_id was not called');

    is_deeply([ $locker->lock_attempts ],
              [ ],
              'Did not try to lock');
    is_deeply([ $locker->unlock_attempts ],
              [ ],
              'Did not try to unlock');
};

subtest 'create with no commit off' => sub {
    plan tests => 5;

    local $ENV{UR_DBI_NO_COMMIT} = 0;
    my $lock_id_called = 0;

    my $locker = LockHijacker->new();

    {    no warnings 'once';
        *TestCreate::lock_id = sub { ++$lock_id_called };
    }
    class TestCreate {
        roles => ['Genome::Role::ObjectWithLockedConstruction'],
    };

    my $o = TestCreate->create();
    ok($o, 'Created object');
    is($lock_id_called, 1, 'method lock_id was called');

    is_deeply([ $locker->lock_attempts ],
        [ { resource_lock => 'TestCreate/1',
            scope => 'site',
            block_sleep => 60,
            wait_announce_interval => 0,
            max_try => 7200,
        } ],
        'One lock attempt during create()');
    is_deeply([ $locker->unlock_attempts ],
              [ ],
              'No unlock attempts after create');

    UR::Context->commit();

    is_deeply([ $locker->unlock_attempts ],
        [ { resource_lock => 'TestCreate/1',
            scope => 'site',
        } ],
        'One unlock attempt during commit');
};

package LockHijacker;

sub new {
    my $class = shift;
    require Genome::Sys::Lock;

    my $self = bless {}, $class;

    $self->{orig_lock_resource} = \&Genome::Sys::Lock::lock_resource;
    $self->{orig_unlock_resource} = \&Genome::Sys::Lock::unlock_resource;
    no warnings 'redefine';
    *Genome::Sys::Lock::lock_resource = \&LockHijacker::_genome_sys_lock_resource;
    *Genome::Sys::Lock::unlock_resource = \&LockHijacker::_genome_sys_unlock_resource;

    return $self;
}

our @lock_attempts;
sub lock_attempts { return @lock_attempts }

sub _genome_sys_lock_resource {
    my($class, %params) = @_;

    push @lock_attempts, \%params;
    return 1;
}

our @unlock_attempts;
sub unlock_attempts { return @unlock_attempts }

sub _genome_sys_unlock_resource {
    my($class, %params) = @_;

    push @unlock_attempts, \%params;
    return 1;
}

sub DESTROY {
    my $self = shift;

    @lock_attempts = @unlock_attempts = ();
    no warnings 'redefine';
    *Genome::Sys::Lock::lock_resource = $self->{orig_lock_resource};
    *Genome::Sys::Lock::unlock_resource = $self->{orig_unlock_resource};
}

