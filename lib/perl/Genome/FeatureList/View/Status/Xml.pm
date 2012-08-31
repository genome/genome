package Genome::FeatureList::View::Status::Xml;

use strict;
use warnings;

use Genome;

class Genome::FeatureList::View::Status::Xml {
    is => 'Genome::View::Status::Xml',
    has_constant => [
        default_aspects => {
            is => 'ARRAY',
            value => [
                'id',
                'name',
                'format',
                'file_content_hash',
                'is_multitracked',
                'is_1_based',
                'source',
                'reference_id',
                {
                    name => 'reference',
                    perspective => 'default',
                    toolkit => 'xml',
                    aspects => [
                        'id', 'data_directory', 'status', 'date_scheduled', 'date_completed', 'name',
                        {
                            name => 'model',
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
                                'user_name',
                            ],
                            subject_class_name => 'Genome::Model',
                        },
                    ],
                },
                'subject_id',
                {
                    name => 'subject',
                    perspective => 'default',
                    toolkit => 'xml',
                    aspects => [
                        'id', 'data_directory', 'status', 'date_scheduled', 'date_completed',
                        {
                            name => 'model',
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
                                'user_name',
                            ],
                            subject_class_name => 'Genome::Model',
                        }
                    ],
                },
                'file_path',
                'content_type',
                {
                    name => 'disk_allocation',
                    perspective => 'default',
                    toolkit => 'xml',
                    aspects => ['absolute_path'],
                },
            ]
        }
    ]
};

1;
