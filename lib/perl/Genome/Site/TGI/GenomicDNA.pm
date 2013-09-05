package Genome::Site::TGI::GenomicDNA;
use strict;
use warnings;
use Genome;

class Genome::Site::TGI::GenomicDNA { 
    id_by => 'gd_id',
    has => [
        'name' => { column_name => 'genomic_name' }, 
    ],
    data_source => 'Genome::DataSource::Oltp',
    table_name => 'genomic_dna',
};


1;

