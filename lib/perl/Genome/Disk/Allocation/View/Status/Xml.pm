package Genome::Disk::Allocation::View::Status::Xml;

use strict;
use warnings;

use Genome;

class Genome::Disk::Allocation::View::Status::Xml {
    is => 'Genome::View::Status::Xml',
    has_constant => [
        default_aspects => {
            is => 'ARRAY',
            value => [
                'absolute_path',
                'kilobytes_requested',
                'owner_class_name',
                'owner_id',
                { name => 'build',
                  perspective => 'default',
                  subject_class_name => 'Genome::Model::Build',
                  toolkit => 'xml',
                  aspects => ['build_id', 'model_id', 'status', 'run_by', 'date_scheduled', 'date_completed' ],
                }
            ],
        }
    ]
};

1;

