package Genome::GoResult;

use strict;
use warnings;
use Genome;

class Genome::GoResult {
    type_name => 'genome go result',
    id_by => [
        go_id => {
           is => 'Number',
        },
    ],
    has => [
        data_directory => {
            is => 'Path',
        },
        chrom_name => {
            is => 'Text',
        },
        interpro_result_id => {
            is => 'Number',
        },
        start => {
            is => 'Number',
        },
        end => {
            is => 'Number',
        },
        transcript_name => { 
            is => 'Number',
        },
        name => {
            is => 'Text',
        },
        term_type => {
            is => 'Text',
        },
    ],
    schema_name => 'files',
    data_source => 'Genome::DataSource::GoResults',
};

1;

