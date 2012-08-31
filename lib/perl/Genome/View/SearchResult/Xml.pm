package Genome::View::SearchResult::Xml;

use strict;
use warnings;
use Genome;

class Genome::View::SearchResult::Xml {
    is => 'UR::Object::View::Default::Xml',
    is_abstract => 1,
    has => [
        solr_doc => {
            is => 'WebService::Solr::Document',
            doc => 'The Solr document that triggered the desire to create this view',
            is_optional => 1,
        }
    ],
    has_constant => [
        perspective => {
            value => 'search-result',
        },
    ],
    doc => 'The base class for creating the XML document representing a "one-line" view of an object'
};

1;
