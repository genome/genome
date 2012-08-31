package Genome::Taxon::View::Status::Xml;

use strict;
use warnings;

use Genome;

class Genome::Taxon::View::Status::Xml {
    is => 'Genome::View::Status::Xml',
    has_constant => [
        default_aspects => {
            is => 'ARRAY',
            value => [
                'id',
                'species_name',
                'species_latin_name',
                'strain_name',
                'estimated_genome_size',
                'domain',
                'gram_stain_category', 
                {
                    name => 'model_member',
                    perspective => 'default',
                    toolkit => 'xml',
                    aspects => [
                        'id',
                        'name',
                        'common_name',
                        {
                            name => 'samples',
                            perspective => 'default',
                            toolkit => 'xml',
                            aspects => [
                                'id',
                                'name',
                                {
                                    name => 'models',
                                    perspective => 'status',
                                    toolkit => 'xml',
                                    subject_class_name => 'Genome::Model',
                                },
                            ]
                        }
                    ]
                },
                {
                    name => 'individuals',
                    perspective => 'default',
                    toolkit => 'xml',
                    aspects => [
                        'id',
                        'name',
                        'common_name',
                        'gender',
# In some cases (e.g. human), including these makes the page take > 30 minutes to generate.
#                        {
#                            name => 'samples',
#                            perspective => 'default',
#                            toolkit => 'xml',
#                            aspects => [
#                                'id',
#                                'name',
#                                {
#                                    name => 'models',
#                                    perspective => 'status',
#                                    toolkit => 'xml',
#                                    subject_class_name => 'Genome::Model',
#                                },
#                            ]
#                        }
                    ]
                },
                {
                    name => 'population_groups',
                    perspective => 'default',
                    toolkit => 'xml',
                    aspects => [
                        'id',
                        'name',
                        'common_name',
# In some cases (e.g. human), including these makes the page take > 30 minutes to generate.
#                        {
#                            name => 'samples',
#                            perspective => 'default',
#                            toolkit => 'xml',
#                            aspects => [
#                                'id',
#                                'name',
#                                {
#                                    name => 'models',
#                                    perspective => 'status',
#                                    toolkit => 'xml',
#                                    subject_class_name => 'Genome::Model',
#                                },
#                            ]
#                        }
                    ]
                },
            ]
        }
    ]
};


1;
