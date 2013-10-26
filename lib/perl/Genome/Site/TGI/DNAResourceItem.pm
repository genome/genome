package Genome::Site::TGI::DNAResourceItem;
use strict;
use warnings;
use Genome;

class Genome::Site::TGI::DNAResourceItem { 
    id_by => 'dri_id',
    has => [
        'name' => { column_name => 'dna_resource_item_name' }, 
    ],
    data_source => 'Genome::DataSource::Oltp',
    table_name => 'dna_resource_item',
};


1;

