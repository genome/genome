package Genome::WorkOrder::View::SearchResult::Xml;

use strict;
use warnings;

use Genome;

class Genome::WorkOrder::View::SearchResult::Xml {
    is => 'Genome::View::SearchResult::Xml',
    has_constant => [
        default_aspects => {
            is => 'ARRAY',
            default => [
                'id',
                'barcode',
                'pipeline',
                'project_name',
                'name',
                'description',
                {
                    name => 'project',
                    aspects => ['id', 'name'],
                    perspective => 'default',
                    toolkit => 'xml'
                },
            ]
        }
    ]
};

1;

