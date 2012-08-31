package Genome::DruggableGene::DrugGeneInteractionReport::Set::View::Status::Xml;

use strict;
use warnings;
use Genome;
use Data::Dumper;
use XML::LibXML;

class Genome::DruggableGene::DrugGeneInteractionReport::Set::View::Status::Xml {
    is => 'Genome::View::Status::Xml',
    has_constant => [
        default_aspects => {
            is => 'ARRAY',
            value => [
                {
                    name => 'members',
                    perspective => 'default',
                    toolkit => 'xml',
                    subject_class_name => 'Genome::DruggableGene::DrugGeneInteractionReport',
                    aspects => [
                        'drug_name',
                        'human_readable_drug_name',
                        'gene_name',
                        'gene_group_name',
                        'interaction_types',
                        'source_db_name',
                        'source_db_version',
                        {
                            name => 'interaction_attributes',
                            perspective => 'default',
                            toolkit => 'xml',
                            aspects => [
                                'name',
                                'value',
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
            ],
        },
    ],
};

1;
