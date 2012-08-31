package Genome::Disk::Volume::View::SearchResult::Xml;

use strict;
use warnings;

use Genome;

class Genome::Disk::Volume::View::SearchResult::Xml {
    is => 'Genome::View::SearchResult::Xml',
    has_constant => [
        default_aspects => {
            is => 'ARRAY',
            default => [
                'dv_id',
                'mount_path',
                'disk_group_names',
            ]
        }
    ]
};

1;
