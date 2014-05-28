package Genome::Config::AnalysisProject::View::Status::Xml;

use strict;
use warnings;

use Genome;

class Genome::Config::AnalysisProject::View::Status::Xml {
    is => 'UR::Object::View::Default::Xml',
    has_constant => [
        default_aspects => {
            is => 'ARRAY',
            value => [
                'id',
                'name',
                'created_by',
                'created_at',
                'updated_at',
            ],
        },
    ],
};

1;
