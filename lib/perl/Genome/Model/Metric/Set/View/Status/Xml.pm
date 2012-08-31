package Genome::Model::Metric::Set::View::Status::Xml;

use strict;
use warnings;

use Genome;

class Genome::Model::Metric::Set::View::Status::Xml {
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
                        'build_id',
                        'name',
                        'value'
                    ],
                    subject_class_name => 'Genome::Model::Metric',
                },
            ]
        }
    ]
};


1;
