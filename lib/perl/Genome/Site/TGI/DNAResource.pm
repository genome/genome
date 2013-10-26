package Genome::Site::TGI::DNAResource;
use strict;
use warnings;
use Genome;

class Genome::Site::TGI::DNAResource { 
    id_by => 'dr_id',
    has => [
        'name' => { column_name => 'dna_resource_prefix' }, 
    ],
    data_source => 'Genome::DataSource::Oltp',
    table_name => 'dna_resource',
};


1;

