package Genome::Disk::Group::View::SearchResult::Xml;

use strict;
use warnings;

use Genome;

class Genome::Disk::Group::View::SearchResult::Xml {
    is => 'Genome::View::SearchResult::Xml',
    has_constant => [
        default_aspects => {
            is => 'ARRAY',
            default => [
                'dg_id',
                'disk_group_name',
                'user_name',
                'group_name',
            ]
        }
    ]
};

1;
