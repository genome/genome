package Genome::Utility::ObjectWithLockedConstruction;

use strict;
use warnings;

use Genome;
use Genome::Sys::LockProxy qw();

class Genome::Utility::ObjectWithLockedConstruction {
    is => 'UR::Object',
    is_abstract => 1,
};

sub create {
    my $class = shift;

    my $obj = $class->get(@_);
    return $obj if $obj;

    if ($ENV{UR_DBI_NO_COMMIT}) {
        return $class->SUPER::create(@_);
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
       $obj = $class->SUPER::create(@_);
       UR::Context->current->add_observer(
           aspect => 'commit',
           callback => sub {
               $lock->unlock();
           }
       );
    }

    return $obj;
}

sub lock_id {
    die("You must implement lock_id in subclasses!");
}
