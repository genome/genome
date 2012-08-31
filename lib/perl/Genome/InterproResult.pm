package Genome::InterproResult;

use strict;
use warnings;
use Genome;

class Genome::InterproResult {
    type_name => 'genome interpro result',
    id_by => [
        interpro_id => {
            is => 'Number',
        },
    ],
    has => [
        chrom_name => {
            is => 'Text',
        },
        transcript_name => { 
            is => 'Text',
        },
        data_directory => {
            is => 'Path',
        },
        start => {
            is => 'Number',
        },
        stop => {
            is => 'Number',
        },
        rid => {
            is => 'Number',
            is_optional => 1,
        },
        setid => {
            is => 'Text',
            is_optional => 1,
        },
        parent_id => {
            is => 'Text',
            is_optional => 1,
        },
        name => {
            is => 'Text',
            is_optional => 1,
        },
        interpro_note => {
            is => 'Text',
            is_optional => 1,
        },
    ],
    schema_name => 'files',
    data_source => 'Genome::DataSource::InterproResults',
};

1;

