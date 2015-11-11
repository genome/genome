package Genome::Role::ObjectWithLockedConstruction;

use strict;
use warnings;

use Genome;
use UR::Role qw(around);
use Genome::Sys::LockProxy qw();

role Genome::Role::ObjectWithLockedConstruction {
    requires => 'lock_id',
};

around 'create' => sub {
    my $orig_create = shift;
    my $class = shift;

    if (my $obj = $class->get(@_)) {
        return $obj;

    } elsif ($ENV{UR_DBI_NO_COMMIT}) {
        return $class->$orig_create(@_);

    } else {
        my $lock_id = $class->lock_id(@_);
        my $lock_var = sprintf('%s/%s', $class, $lock_id);

        return $class->_create_with_lock($lock_var, $orig_create, @_);
    }
};

sub _create_with_lock {
    my $class = shift;
    my $lock_var = shift;
    my $orig_create = shift;

    my $lock = Genome::Sys::LockProxy->new(
        resource => $lock_var,
        scope => 'site',
    )->lock();
    die("Unable to get lock!") unless $lock;

    my $obj = $class->load(@_);
    unless ($obj) {
        $obj = $class->$orig_create(@_);
        Genome::Sys::CommitAction->create(
            on_commit => sub {
                $lock->unlock();
            }
        );
    }

    return $obj;
}

1;
