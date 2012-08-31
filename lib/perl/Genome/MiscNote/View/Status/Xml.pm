package Genome::MiscNote::View::Status::Xml;

use strict;
use warnings;

use Genome;

class Genome::MiscNote::View::Status::Xml {
    is => 'Genome::View::Status::Xml',
    has_constant => [
        default_aspects => {
            is => 'ARRAY',
            value => [
                'editor_id',
                'entry_date',
                'entry_date_sort',
                'header_text',
                'body_text',
                {
                    name => 'subject',
                    perspective => 'default',
                    toolkit => 'xml'
                },
            ],
        }
    ]
};
1;
