package Genome::View::SearchResult::Html;

use strict;
use warnings;

use Genome;

class Genome::View::SearchResult::Html {
    is => 'UR::Object::View::Default::Html',
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
    doc => 'Concrete classes that want to be able to transform their SearchResult::Xml views to html should inherit from this class.'
};

1;
