package Genome::Task::Set::View::Status::Xml;

use strict;
use warnings;

use Genome;

class Genome::Task::Set::View::Status::Xml {
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
                        'command_class',
                        'time_submitted',
                        'time_started',
                        'time_finished',
                        'status',
                        'user_id',
                    ],
                    subject_class_name => 'Genome::Task',
                },
            ]
        }
    ]
};


1;
