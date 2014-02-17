package Genome::Model::Build::Set::View::Status::Xml;

use strict;
use warnings;

use Genome;

class Genome::Model::Build::Set::View::Status::Xml {
    is => 'Genome::View::Status::Xml',
    has_constant => [
        default_aspects => {
            is => 'ARRAY',
            value => [
                'rule_display',
                {
                    name => 'members',
                    perspective => 'default',
                    toolkit => 'xml',
                    aspects => [
                        'id',
                        'data_directory',
                        'status',
                        'date_scheduled',
                        'date_completed',
                        'run_by',
                        {
                            name => 'model',
                            aspects => [
                                'genome_model_id',
                                'name',
                                {
                                    name => 'subject',
                                    aspects => ['id','name','subclass_name'],
                                    perspective => 'default',
                                    toolkit => 'xml',
                                    subject_class_name => 'Genome::Subject',
                                },
                                'creation_date',
                                'created_by',
                                'run_as',
                            ],
                            perspective => 'default',
                            toolkit => 'xml',
                            subject_class_name => 'Genome::Model',
                        }
                    ],
                    subject_class_name => 'Genome::Model::Build',
                },
            ]
        }
    ]
};


1;
