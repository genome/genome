package Genome::Utility::ObjectWithTimestamps;

use strict;
use warnings;

use Genome;

class Genome::Utility::ObjectWithTimestamps {
    subclass_description_preprocessor => __PACKAGE__ . '::_preprocess_subclass_description',
    is => 'UR::Object',
    is_abstract => 1,
};

Genome::Utility::ObjectWithTimestamps->add_observer(callback => \&is_updated);
Genome::Utility::ObjectWithTimestamps->add_observer(aspect => 'create',  callback => \&is_created);

sub is_created {
    my $self = shift;
    unless ($self->created_at) {
        $self->created_at(UR::Context->current->now);
    }
    return 1;
}

sub is_updated {
    my ($self, $aspect) = @_;
    if (ref($self) && $aspect ne 'commit' && $aspect ne 'load') {
        $self->updated_at(UR::Context->current->now);
    }
}

sub _preprocess_subclass_description {
    my ($class, $desc) = @_;
    for ('created_at', 'updated_at') {
        $desc->{has}{$_} = _property_hash_for_name($_);
    }
    return $desc;
}

sub _property_hash_for_name {
    return {
        property_name => shift,
        is => 'Timestamp'
    };
}

1;