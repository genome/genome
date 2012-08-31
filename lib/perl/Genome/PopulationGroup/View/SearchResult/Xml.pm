package Genome::PopulationGroup::View::SearchResult::Xml;

use strict;
use warnings;

use Genome;

class Genome::PopulationGroup::View::SearchResult::Xml {
    is => 'Genome::View::SearchResult::Xml',
    has_constant => [
        default_aspects => {
            is => 'ARRAY',
            default => [
                'id',
                'name',
                'common_name',
                'description',
                {
                    name => 'members',
                    aspects => ['id', 'name', 'common_name', ],
                    perspective => 'default',
                    toolkit => 'xml'
                },
            ]
        }
    ]
};

1;
