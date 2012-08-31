package Genome::Sample::View::Detail::Xml;

use strict;
use warnings;

use Genome;

class Genome::Sample::View::Detail::Xml {
    is => 'UR::Object::View::Default::Xml',
    has_constant => [
        default_aspects => {
            is => 'ARRAY',
            value => [
                'id',
                'name',
                {
                    name => 'attributes',
                    perspective => 'default',
                    toolkit => 'xml',
                    aspects => ['nomenclature_name','nomenclature_field_name','attribute_label','attribute_value','subject_id']
                }
            ]
        }
    ]
};


1;
