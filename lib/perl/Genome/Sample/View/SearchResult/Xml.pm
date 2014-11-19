package Genome::Sample::View::SearchResult::Xml;

use strict;
use warnings;

use Genome;

class Genome::Sample::View::SearchResult::Xml {
    is => 'Genome::View::SearchResult::Xml',
    has_constant => [
        default_aspects => {
            is => 'ARRAY',
            default => [
                'id',
                'name',
                'common_name',
                'extraction_label',
                'extraction_type',
                'extraction_desc',
                'tissue_label',
                'tissue_desc',
                'organ_name',
                {
                    name => 'taxon',
                    aspects => ['id', 'species_name'],
                    perspective => 'default',
                    toolkit => 'xml'
                },
                {
                    name => 'source',
                    aspects => ['id', 'name', 'common_name'],
                    perspective => 'default',
                    toolkit => 'xml'
                }
            ]
        }
    ]
};

1;
