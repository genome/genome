package Genome::Site::TGI::IPRProduct;
use strict;
use warnings;
use Genome;

class Genome::Site::TGI::IPRProduct { 
    table_name => '(select * from gsc.dna d join gsc.ipr_product i on i.ipr_id = d.dna_id) d',
    id_by => 'dna_id',
    has => [
        'name' => { column_name => 'dna_name' }, 
    ],
    data_source => 'Genome::DataSource::Oltp'
};


1;

