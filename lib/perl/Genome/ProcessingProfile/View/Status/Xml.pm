package Genome::ProcessingProfile::View::Status::Xml;

use strict;
use warnings;

use Genome;

class Genome::ProcessingProfile::View::Status::Xml {
    is => 'Genome::View::Status::Xml',
    has_constant => [
        default_aspects => {
            is => 'ARRAY',
            value => [
                'id',
                'name',
                'type_name',
                'supersedes',
                {
                    name => 'params',
                    perspective => 'default',
                    toolkit => 'xml',
                    aspects => [
                        'id', 'name', '_value_scalar_or_object',
                    ]
                }
            ]
        }
    ]
};


1;
