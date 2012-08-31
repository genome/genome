package Genome::Site::TGI::PCRProduct;
use strict;
use warnings;
use Genome;

class Genome::Site::TGI::PCRProduct { 
    table_name => 'GSC.PCR_PRODUCT@oltp pcr_product',
    id_by => 'pcr_id',
    has => [
        'name' => { column_name => 'pcr_name' }, 
    ],
    data_source => 'Genome::DataSource::GMSchema'
};


1;

