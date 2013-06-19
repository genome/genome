package Genome::Utility::ObjectWithTimestamps;

use strict;
use warnings;

use Genome;

#If the object mixing this in has an underlying database table, you will need to
#explicitly define created_at and updated_at properties in your class in order
#to have your created_at and updated_at fields synced to the the database
class Genome::Utility::ObjectWithTimestamps {
    is => 'UR::Object',
    has_optional => [
        created_at => {
            is => 'Timestamp',
        },
        updated_at => {
            is => 'Timestamp',
        },
    ],
};

Genome::Utility::ObjectWithTimestamps->add_observer(callback => \&is_updated);

sub create {
    my $class = shift;
    my $self = $class->SUPER::create(@_);
    eval {
        $self->created_at(UR::Context->current->now);
    };
    if(my $error = $@) {
        $self->delete();
        die($error);
    }
    return $self;
}

sub is_updated {
    my ($self, $aspect) = @_;
    if (ref($self) && $aspect ne 'commit' && $aspect ne 'load') {
        $self->updated_at(UR::Context->current->now);
    }
}

1;