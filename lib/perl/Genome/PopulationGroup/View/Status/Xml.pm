package Genome::PopulationGroup::View::Status::Xml;

use strict;
use warnings;

use Genome;

class Genome::PopulationGroup::View::Status::Xml {
    is => 'Genome::View::Status::Xml',
    has_constant => [
        default_aspects => {
            is => 'ARRAY',
            value => [
                'id',
                'name',
                {
                    name => 'members',
                    perspective => 'default',
                    toolkit => 'xml',
                    aspects => [
                        'id',
                        'name',
                        'common_name',
                        'gender',
                        {
                            name => 'samples',
                            perspective => 'default',
                            toolkit => 'xml',
                            aspects => [
                                'id',
                                'name',
                                {
                                    name => 'models',
                                    perspective => 'default',
                                    toolkit => 'xml',
                                    aspects => [
                                        'genome_model_id',
                                        'name',
                                        {
                                            name => 'subject',
                                            aspects => ['id', 'name', 'subclass_name'],
                                            perspective => 'default',
                                            toolkit => 'xml',
                                        },
                                        {
                                            name => 'processing_profile',
                                            aspects => ['id', 'name'],
                                            perspective => 'default',
                                            toolkit => 'xml'
                                        },
                                        'creation_date',
                                        'created_by',
                                        'run_as',
                                        {
                                            name => 'last_succeeded_build',
                                            aspects => [ 'id', 'data_directory', 'status', ],
                                            perspective => 'default',
                                            toolkit => 'xml',
                                            subject_class_name => 'Genome::Model::Build',
                                        }
                                    ],
                                    subject_class_name => 'Genome::Model',
                                }
                            ]
                        }
                    ]  
                },
            ]
        }
    ]
};


1;
