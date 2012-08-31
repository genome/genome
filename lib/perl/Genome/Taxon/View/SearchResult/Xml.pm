package Genome::Taxon::View::SearchResult::Xml;

use strict;
use warnings;

use Genome;

class Genome::Taxon::View::SearchResult::Xml {
    is => 'Genome::View::SearchResult::Xml',
    has_constant => [
        default_aspects => {
            is => 'ARRAY',
            default => [
                'id',
                'domain',
                'species_name',
                'domain',
                'species_latin_name',
                'strain_name',
                'ncbi_taxon_id',
            ]
        }
    ]
};

1;
