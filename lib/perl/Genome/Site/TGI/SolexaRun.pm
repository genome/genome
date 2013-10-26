package Genome::Site::TGI::SolexaRun;
use strict;
use warnings;
use Genome;

class Genome::Site::TGI::SolexaRun { 
    id_by => 'er_id',
    has => [
        'name' => { column_name => 'flow_cell_id' }, 
    ],
    data_source => 'Genome::DataSource::Oltp',
    table_name => 'solexa_run',
};


1;

