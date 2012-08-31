package Genome::WorkOrder::View::Status::Xml;

use strict;
use warnings;

use Genome;

class Genome::WorkOrder::View::Status::Xml {
    is => 'Genome::View::Status::Xml',
    has_constant => [
        default_aspects => {
            is => 'ARRAY',
            value => [
                'id',
                'pipeline',
                {
                    name => 'project',
                    perspective => 'default',
                    toolkit => 'xml',
                    subject_class_name => 'Genome::Site::TGI::Project',
                    aspects => ['name', 'status', 'description', 'project_type'],
                },
                {
                    name => 'items',
                    perspective => 'Status',
                    toolkit => 'xml',
                    subject_class_name => 'Genome::WorkOrderItem',
                }
            ],
        }
    ]
};



1;
