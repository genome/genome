package Genome::DruggableGene::DrugNameReport::Set::View::Status::Xml;

use strict;
use warnings;
use Genome;
use Data::Dumper;
use XML::LibXML;

class Genome::DruggableGene::DrugNameReport::Set::View::Status::Xml {
    is => 'Genome::View::Status::Xml',
    has_constant => [
        default_aspects => {
            is => 'ARRAY',
            value => [
                {
                    name => 'members',
                    perspective => 'default',
                    toolkit => 'xml',
                    subject_class_name => 'Genome::DruggableGene::DrugNameReport',
                    aspects => [
                        'id',
                        'name',
                        'source_db_name',
                        'source_db_version',
                        'source_db_url',
                        'original_data_source_url',
                        'human_readable_name',
                        {
                            name => 'drug_alt_names',
                            perspective => 'default',
                            toolkit => 'xml',
                            aspects => [
                                'alternate_name',
                                'nomenclature',
                            ],
                        },
                        {
                            name => 'drug_categories',
                            perspective => 'default',
                            toolkit => 'xml',
                            aspects => [
                                'category_name',
                                'category_value',
                            ],
                        },
                        {
                            name => 'citation',
                            perspective => 'default',
                            toolkit => 'xml',
                            aspects => [
                                'citation',
                            ],
                        },
                    ],
                },
                'name',
            ]
        }
    ],
};

1;
