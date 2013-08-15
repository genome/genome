package Genome::Site::TGI::PCRProduct;
use strict;
use warnings;
use Genome;

class Genome::Site::TGI::PCRProduct { 
    id_by => 'pcr_id',
    has => [
        'name' => { column_name => 'pcr_name' }, 
    ],
    data_source => 'Genome::DataSource::Oltp',
    table_name => 'pcr_product',
};


1;

