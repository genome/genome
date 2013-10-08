package Genome::Timeline::Event;

use strict;
use warnings;

use Genome;
use Sub::Install ();

class Genome::Timeline::Event {
    is => ['Genome::Utility::ObjectWithTimestamps', 'Genome::Utility::ObjectWithCreatedBy'],
    is_abstract => 1,
    id_generator => '-uuid',
    data_source => 'Genome::DataSource::GMSchema',
    id_by => [
        id => {
            is => 'Text',
            len => 64,
        },
    ],
    has => [
        object => {
            is => 'UR::Object',
            id_by => 'object_id',
            id_class_by => 'object_class_name',
        },
        name => {
            is => 'Text',
        },
        object_id => {
            is => 'Text',
        },
        object_class_name => {
            is => 'Text',
        },
        reason => {
            is => 'Text',
        },
    ],
};

sub _add {
    my $class = shift;
    my ($name, $reason, $instance) =  @_;

    my %prop_hash = $class->_properties_to_snapshot();

    return $class->create(
        name => $name,
        object => $instance,
        reason => $reason,
        map { $prop_hash{$_} => $instance->$_ } keys %prop_hash,
    );
}

#hash mapping value names on the object (key) to database columns for the event (value)
sub _properties_to_snapshot { die('Please provide a map of properties to be saved!'); }

sub _define_event_constructors {
    my ($class, $into, $values) = @_;
    for my $event_name (@$values) {
        Sub::Install::install_sub({
                into => $into,
                as => $event_name,
                code => sub { shift->_add($event_name, @_); }
        });
    }
}

1;
