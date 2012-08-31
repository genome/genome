package Genome::Wiki::Document::View::SearchResult::Xml;

use strict;
use warnings;

use Genome;

class Genome::Wiki::Document::View::SearchResult::Xml {
    is => 'Genome::View::SearchResult::Xml',
    has_constant => [
        default_aspects => {
            is => 'ARRAY',
            default => [
                'id',
                'title',
                'user',
                'timestamp',
                'content',
            ]
        }
    ]
};




1;
