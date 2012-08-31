package Genome::Sys::User::View::Status::Xml;

use strict;
use warnings;

use Genome;

class Genome::Sys::User::View::Status::Xml {
    is => 'Genome::View::Status::Xml',
    has_constant => [
        default_aspects => {
            is => 'ARRAY',
            value => [
                'id',
                'name',
                'username',
                'email',
                {
                    name => 'projects',
                    toolkit => 'xml',
                    perspective => 'default',
                    aspects => [ 'id', 'name', 'parts_count' ]
                }
            ]
        }
    ]
};


1;
