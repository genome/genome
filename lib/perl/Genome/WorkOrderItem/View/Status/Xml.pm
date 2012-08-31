package Genome::WorkOrderItem::View::Status::Xml;

use strict;
use warnings;

use Genome;

class Genome::WorkOrderItem::View::Status::Xml {
    is => 'Genome::View::Status::Xml',
    has_constant => [
        default_aspects => {
            is => 'ARRAY',
            value => [
                'id',
                {
                    name => 'event_statuses',
                    subject_class_name => 'UR::Value::HASH',
                    perspective => 'default',
                    toolkit => 'xml',
                },
                {
                    name => 'sample',
                    perspective => 'default',
                    toolkit => 'xml',
                    subject_class_name => 'Genome::Sample',
                    aspects => [
                        'name',
                        'common_name',
                        'species_name'
                    ]
                },
                {
                    name => 'models',
                    perspective => 'default',
                    toolkit => 'xml',
                    subject_class_name => 'Genome::Model',
                    aspects => [
                        'id',
                        'name',
                        {
                            name => 'latest_build',
                            perspective => 'default',
                            toolkit => 'xml',
                            subject_class_name => 'Genome::Model::Build',
                            aspects => [
                                'id',
                                'master_event_status'
                            ]
                        }
                    ]
                }
            ]
        }
    ]
};



1;
