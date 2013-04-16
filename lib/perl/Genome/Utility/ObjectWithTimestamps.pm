package Genome::Utility::ObjectWithTimestamps;

use strict;
use warnings;

use Genome;

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
    if (ref($self) && $aspect ne 'commit') {
        $self->updated_at(UR::Context->current->now);
    }
}

1;