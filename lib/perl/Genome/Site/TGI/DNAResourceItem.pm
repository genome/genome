package Genome::Site::TGI::DNAResourceItem;
use strict;
use warnings;
use Genome;

class Genome::Site::TGI::DNAResourceItem { 
    table_name => 'GSC.DNA_RESOURCE_ITEM@oltp dna_resource_item',
    id_by => 'dri_id',
    has => [
        'name' => { column_name => 'dna_resource_item_name' }, 
    ],
    data_source => 'Genome::DataSource::GMSchema'
};


1;

