package Genome::Disk::Volume::View::Status::Xml;

use strict;
use warnings;

use Genome;

class Genome::Disk::Volume::View::Status::Xml {
    is => 'Genome::View::Status::Xml',
    has_constant => [
        default_aspects => {
            is => 'ARRAY',
            value => [
                'dv_id',
                'mount_path',
                'disk_status',
                'can_allocate',
                'unallocated_kb',
                'total_kb',
                'disk_group_names',
                {        name => 'allocations',
                  perspective => 'status',
                      toolkit => 'xml',
                }
            ],
        }
    ]
};

1;

