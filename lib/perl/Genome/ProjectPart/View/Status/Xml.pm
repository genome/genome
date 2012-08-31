package Genome::ProjectPart::View::Status::Xml;

use strict;
use warnings;

use Genome;

class Genome::ProjectPart::View::Status::Xml {
    is => 'Genome::View::Status::Xml',
    has => [
        default_aspects => {
            is => 'ARRAY',
            value => [
                'entity_id',
                'entity_class_name',
                'entity_class_name_pretty',
                'entity',
            ]
        }
    ]
};


1;
