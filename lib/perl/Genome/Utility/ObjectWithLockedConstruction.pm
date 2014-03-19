package Genome::Utility::ObjectWithLockedConstruction;

use strict;
use warnings;

use Genome;

class Genome::Utility::ObjectWithLockedConstruction {
    is => 'UR::Object',
    is_abstract => 1,
};

sub create {
    my $class = shift;
    if ($ENV{UR_DBI_NO_COMMIT}) {
        return $class->SUPER::create(@_);
    } else {
        my $lock_id = $class->lock_id(@_);

        my $lock_var = sprintf('%s/%s/%s',
            $ENV{GENOME_LOCK_DIR},
            $class,
            $lock_id);

        my $obj = $class->get(@_);
        if($obj) {
            return $obj;
        } else {
            return $class->create_with_lock($lock_var, @_);
        }
    }
}

sub create_with_lock {
    my $class = shift;
    my $lock_var = shift;

    my $lock = Genome::Sys->lock_resource(resource_lock => $lock_var);
    die("Unable to get lock!") unless $lock;

    my $obj = $class->load(@_);
    if($obj) {
        return $obj;
    } else {
       $obj = $class->SUPER::create(@_);
       UR::Context->current->add_observer(
           aspect => 'commit',
           callback => sub {
               Genome::Sys->unlock_resource(resource_lock => $lock_var);
           }
       );
    }

    return $obj;
}

sub lock_id {
    die("You must implement lock_id in subclasses!");
}