package Genome::Individual::View::SearchResult::Xml;

use strict;
use warnings;

use Genome;

class Genome::Individual::View::SearchResult::Xml {
    is => 'Genome::View::SearchResult::Xml',
    has_constant => [
        default_aspects => {
            is => 'ARRAY',
            default => [
                'common_name',
                'name',
                'upn',
                'gender',
                'description',
                'species_name',
            ]
        }
    ]
};

1;
