package Genome::Model::View::SearchResult::Xml;

use strict;
use warnings;

use Genome;

class Genome::Model::View::SearchResult::Xml {
    is => 'Genome::View::SearchResult::Xml',
    has_constant => [
        default_aspects => {
            is => 'ARRAY',
            default => [
                'genome_model_id',
                'name',
                'subject_id',
                'subject_class_name',
                'data_directory',
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
                    name => 'last_complete_build',
                    aspects => [
                        'id', 'data_directory'
                    ],
                    perspective => 'default',
                    toolkit => 'xml',
                    subject_class_name => 'Genome::Model::Build',
                }
            ]
        }
    ]
};

1;
