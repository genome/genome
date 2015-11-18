package Genome::Timeline::Event::Base;

use strict;
use warnings;

use Genome;

class Genome::Timeline::Event::Base {
    is => 'Genome::Timeline::Event',
    roles => ['Genome::Role::ObjectWithCreatedBy', 'Genome::Role::ObjectWithTimestamps'],
    table_name => 'timeline.base',
};


1;
