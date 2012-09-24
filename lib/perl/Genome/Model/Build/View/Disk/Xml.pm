package Genome::Model::Build::View::Disk::Xml;

use strict;
use warnings;

use Genome;

class Genome::Model::Build::View::Disk::Xml {
    is => 'Genome::View::Status::Xml',
    has_constant => [
        default_aspects => {
            is => 'ARRAY',
            value => [
                'id',
                'model_id'
            ]
        }
    ]
};


1;
