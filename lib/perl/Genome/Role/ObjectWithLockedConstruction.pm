package Genome::Role::ObjectWithLockedConstruction;

use strict;
use warnings;

use UR::Role;
use Genome::Sys::LockProxy qw();

role Genome::Role::ObjectWithLockedConstruction {
    requires => 'lock_id',
};

my %super_create;
sub create {
    my $class = shift;

    my $obj = $class->get(@_);
    return $obj if $obj;

    if ($ENV{UR_DBI_NO_COMMIT}) {
        my $super_create = $super_create{$class} ||= $class->super_can('create');
        return $class->$super_create(@_);
    } else {
        my $lock_id = $class->lock_id(@_);
        my $lock_var = sprintf('%s/%s', $class, $lock_id);

        return $class->create_with_lock($lock_var, @_);
    }
}

sub create_with_lock {
    my $class = shift;
    my $lock_var = shift;

    my $lock = Genome::Sys::LockProxy->new(
        resource => $lock_var,
        scope => 'site',
    )->lock();
    die("Unable to get lock!") unless $lock;

    my $obj = $class->load(@_);
    if($obj) {
        return $obj;
    } else {
        my $super_create = $super_create{$class} ||= $class->super_can('create');
        $obj = $class->$super_create(@_);
        Genome::Sys::CommitAction->create(
            on_commit => sub {
                $lock->unlock();
            }
        );
    }

    return $obj;
}

1;
