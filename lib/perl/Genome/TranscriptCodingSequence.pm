package Genome::TranscriptCodingSequence;

use strict;
use warnings;
use Genome;

class Genome::TranscriptCodingSequence {
    type_name => 'genome transcript coding sequence',
    id_by => [
        transcript_id => {
            is => 'Text',
        },
    ],
    has => [
        sequence => {
            is => 'Text',
        },
        data_directory => {
            is => 'Path',
        },
    ],
    schema_name => 'files',
    data_source => 'Genome::DataSource::TranscriptCodingSequences',
};

1;

