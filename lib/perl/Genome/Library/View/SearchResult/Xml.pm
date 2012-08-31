package Genome::Library::View::SearchResult::Xml;

use strict;
use warnings;

use Genome;

class Genome::Library::View::SearchResult::Xml {
    is => 'Genome::View::SearchResult::Xml',
    has_constant => [
        default_aspects => {
            is => 'ARRAY',
            default => [
                'library_id',
                'name',
                {
                    name => 'sample',
                    aspects => [ 'id', 'name', ],
                    perspective => 'default',
                    toolkit => 'xml'
                },
                {
                    name => 'taxon',
                    aspects => [ 'id', 'species_name', ],
                    perspective => 'default',
                    toolkit => 'xml'
                },
            ]
        }
    ]
};

1;
