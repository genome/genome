package Genome::SoftwareResult::View::Status::Xml;

use strict;
use warnings;

use Genome;

class Genome::SoftwareResult::View::Status::Xml {
    is => 'Genome::View::Status::Xml',
    has_constant => [
        default_aspects => {
            is => 'ARRAY',
            value => [
                'id',
                'subclass_name',
                'test_name',
                'output_dir',
                {
                    name => 'params',
                    perspective => 'default',
                    toolkit => 'xml',
                    aspects => [
                        'name',
                        'value_id',
                    ],
                },
                {
                    name => 'inputs',
                    perspective => 'default',
                    toolkit => 'xml',
                    aspects => [
                        'name',
                        'value_id',
                    ],
                },
                {
                    name => 'metrics',
                    perspective => 'default',
                    toolkit => 'xml',
                    aspects => [
                        'metric_name',
                        'metric_value',
                    ],
                },
                {
                    name => 'users',
                    perspective => 'default',
                    toolkit => 'xml',
                    aspects => [
                        'label',
                        'user_id',
                        'user_class',
                        {
                            name => 'user',
                            perspective => 'default',
                            toolkit => 'xml',
                            aspects => [
                                'id',
                                '__display_name__',
                            ],
                        },
                    ],
                },

            ]
        },
    ]
};


1;
