package Genome::Site::TGI::GenomicDNA;
use strict;
use warnings;
use Genome;

class Genome::Site::TGI::GenomicDNA { 
    table_name => 'GSC.GENOMIC_DNA@oltp genomic_dna',
    id_by => 'gd_id',
    has => [
        'name' => { column_name => 'genomic_name' }, 
    ],
    data_source => 'Genome::DataSource::GMSchema'
};


1;

