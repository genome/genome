package Genome::Project::View::Status::Xml;

use strict;
use warnings;

use Genome;

class Genome::Project::View::Status::Xml {
    is => 'Genome::View::Status::Xml',
    has_constant => [
        default_aspects => {
            is => 'ARRAY',
            value => [
                'id',
                'name',
                {
                    name => 'parts',
                    toolkit => 'xml',
                    perspective => 'status',
                }
            ]
        }
    ]
};


1;
