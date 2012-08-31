package Genome::Nomenclature::View::Detail::Xml;

use strict;
use warnings;

use Genome;

class Genome::Nomenclature::View::Detail::Xml {
    is => 'UR::Object::View::Default::Xml',
    has_constant => [
        default_aspects => {
            is => 'ARRAY',
            value => [
                'id',
                'name',
                { name => 'fields',
                    perspective => 'default',
                    toolkit => 'xml', 
                    aspects => ['id', 'name', 
                        { name => 'enumerated_values',
                            perspective => 'default',
                            toolkit => 'xml', 
                            aspects => ['id', 'value']
                        }
                    ],
                }
            ]
        }
    ]
};


1;
