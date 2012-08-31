package Genome::Library::View::Status::Xml;

use strict;
use warnings;

use Genome;

class Genome::Library::View::Status::Xml {
    is => 'Genome::View::Status::Xml',
    has_constant => [
        default_aspects => {
            is => 'ARRAY',
            value => [
                'id',
                'name',
                {
                    name => 'models',
                    perspective => 'status',
                    toolkit => 'xml',
                    subject_class_name => 'Genome::Model',
                },
                {
                    name => 'sample',
                    perspective => 'default',
                    toolkit => 'xml',
                    aspects => [
                        'id',
                        'name',
                    ]
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
                        'ncbi_taxon_id',
                    ]
                }
            ]
        }
    ]
};


1;
