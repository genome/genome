package Genome::InstrumentData::FlowCell::View::SearchResult::Xml;

use strict;
use warnings;

use Genome;

class Genome::InstrumentData::FlowCell::View::SearchResult::Xml {
    is => 'Genome::View::SearchResult::Xml',
    has_constant => [
        default_aspects => {
            is => 'ARRAY',
            default => [
                'flow_cell_id',
                'machine_name',
                'run_name',
                'run_type',
                {
                    name => 'lanes',
                    perspective => 'default',
                    toolkit => 'xml',
                    aspects => [
                        'id',
                    ],
                },
            ]
        }
    ]
};

1;
