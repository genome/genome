package Genome::ProcessingProfile::View::SearchResult::Xml;

use strict;
use warnings;

use Genome;

class Genome::ProcessingProfile::View::SearchResult::Xml {
    is => 'Genome::View::SearchResult::Xml',
    has_constant => [
        default_aspects => {
            is => 'ARRAY',
            default => [
                'id',
                'name',
                'type_name',
            ]
        }
    ]
};

1;
