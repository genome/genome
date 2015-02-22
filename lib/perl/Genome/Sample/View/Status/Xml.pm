package Genome::Sample::View::Status::Xml;

use strict;
use warnings;

use Genome;

class Genome::Sample::View::Status::Xml {
    is => 'Genome::View::Status::Xml',
    has_constant => [
        default_aspects => {
            is => 'ARRAY',
            value => [
                'id',
                'name',
                'common_name',
                'extraction_label',
                'extraction_type',
                'extraction_desc',
                'tissue_label',
                'tissue_desc',
                'organ_name',
                'individual_common_name',
                {
                    name => 'models',
                    perspective => 'status',
                    toolkit => 'xml',
                    subject_class_name => 'Genome::Model',
                },
                {
                    name => 'taxon',
                    perspective => 'default',
                    toolkit => 'xml',
                    aspects => [
                        'id',
                        'species_name',
                        'species_latin_name',
                        'strain_name',
                    ]
                },
                {
                    name => 'source',
                    perspective => 'default',
                    toolkit => 'xml',
                    aspects => [
                        'id',
                        'name',
                        'common_name',
                    ]
                },
                {
                    name => 'libraries',
                    perspective => 'default',
                    toolkit => 'xml',
                    aspects => [
                        'id',
                        'name',
                    ]
                },
            ]
        }
    ]
};


1;
