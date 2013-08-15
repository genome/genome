package Genome::Site::TGI::Ligation;
use strict;
use warnings;
use Genome;

class Genome::Site::TGI::Ligation { 
    id_by => 'lig_id',
    has => [
        'name' => { column_name => 'ligation_name' }, 
    ],
    data_source => 'Genome::DataSource::Oltp',
    table_name => 'ligations',
};


1;

