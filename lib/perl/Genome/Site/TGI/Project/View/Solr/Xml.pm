package Genome::Site::TGI::Project::View::Solr::Xml;

use strict;
use warnings;

use Genome;

class Genome::Site::TGI::Project::View::Solr::Xml {
    is => 'Genome::View::Solr::Xml',
    has_field => [
        type => {
            is => 'Text',
            default => 'project'
        },
        default_aspects => {
            is => 'ARRAY',
            default => [
                {
                    name => 'description',
                    position => 'content',
                },
                {
                    name => 'status',
                    position => 'content',
                }
            ],
        }
    ]
};

1;
