package Genome::Timeline::Event::AllocationEventType;

use strict;
use warnings;

use Genome;

class Genome::Timeline::Event::AllocationEventType {
    is => 'UR::Object',
    table_name => 'timeline.allocation_event_type',
    data_source => 'Genome::DataSource::GMSchema',
    id_by => [
        id => {
            is => 'Text',
            valid_values => [
               'created',
               'purged',
               'preserved',
               'moved',
               'reallocated',
               'archived',
               'strengthened',
               'weakened',
               'unarchived',
               'unpreserved',
               'finalized',
               'invalidated',
            ]
        }
    ],
    has => [
        events => {
            is_many => 1,
            is => 'Genome::Timeline::Event::Allocation',
            reverse_as => 'type',
        }
    ],
};

1;
