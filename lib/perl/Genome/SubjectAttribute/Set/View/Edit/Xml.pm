package Genome::SubjectAttribute::Set::View::Edit::Xml;

use strict;
use warnings;

use Genome;

class Genome::SubjectAttribute::Set::View::Edit::Xml {
    is => 'UR::Object::Set::View::Default::Xml',
    has_constant => [
        default_aspects => {
            is => 'ARRAY',
            value => [
                'id',
                {
                    name => 'members',
                    perspective => 'default',
                    toolkit => 'xml',
                    subject_class_name => 'Genome::SubjectAttribute',
                    aspects => [
                        'id',
                        'subject_name',
                        'nomenclature_field_name',
                        'attribute_value',
                        {
                            name => 'all_nomenclature_fields',
                            perspective => 'default',
                            toolkit => 'xml',
                            subject_class_name => 'Genome::Nomenclature::Field',
                            aspects => ['name','type','nomenclature_name',
                                {   name => 'enumerated_values',
                                    subject_class_name => 'Genome::Nomenclature::Field::EnumValue',
                                    perspective => 'default',
                                    toolkit => 'xml',
                                    aspects => ['value']
                            }]
                        }
                    ]
                }
            ]
        }
    ]
};


1;


