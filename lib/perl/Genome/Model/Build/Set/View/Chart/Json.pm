package Genome::Model::Build::Set::View::Chart::Json;

use strict;
use warnings;

use Genome;

class Genome::Model::Build::Set::View::Chart::Json {
    is => 'UR::Object::Set::View::Default::Json',
    has_constant => [
        default_aspects => {
            is => 'ARRAY',
            value => [
                'rule_display',
                {
                    name => 'members',
                    subject_class_name => 'Genome::Model::Build',
                    perspective => 'default',
                    toolkit => 'json',
                    aspects => [
                        {
                          name => 'metrics',
                          subject_class_name => 'Genome::Model::Metric',
                          perspective => 'default',
                          toolkit => 'json',
                          aspects => [
                            'build_id',
                            'name',
                            'value'
                          ],
                        }
                    ],
                },
            ]
        }
    ]
};

1;
