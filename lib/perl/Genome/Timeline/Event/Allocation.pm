package Genome::Timeline::Event::Allocation;

use strict;
use warnings;

use Genome;

class Genome::Timeline::Event::Allocation {
    is => 'Genome::Timeline::Event',
    table_name => 'timeline.allocation',
    id_by => [
        id => {
            is => 'Text',
            len => 64,
        },
    ],
    has => [
        object_class_name => {
            is_constant => 1,
            default_value => 'Genome::Disk::Allocation',
            valid_values => ['Genome::Disk::Allocation'],
        },
        #allocation => {
            #via => '__self__',
            #to => 'object',
        #},
        kilobytes_requested => {
            is => 'Number',
        },
        absolute_path => {
            is => 'Text',
        },
        status => {
            is => 'Text',
            valid_values => Genome::Disk::Allocation->__meta__->property('status')->valid_values()
        },
        type => {
            is => 'Genome::Timeline::Event::AllocationEventType',
            id_by => 'name',
        },
    ],
};

Genome::Timeline::Event->_define_event_constructors(
    __PACKAGE__,
    Genome::Timeline::Event::AllocationEventType->__meta__->property('id')->valid_values(),
);

sub _properties_to_snapshot {
    return (
        status => 'status',
        absolute_path => 'absolute_path',
        kilobytes_requested => 'kilobytes_requested',
    );
}

1;
