package Genome::Site::TGI::Gene; 

use strict;
use warnings;

use Genome;

class Genome::Site::TGI::Gene {
    table_name => 'GSC.GENE',
    id_by => [
        gene_id => { 
            is => 'Number', 
            column_name => 'GENE_ID', 
        },
    ],
    has_optional => [
        gene_name => {
            is => 'Text', 
            column_name => 'GENE_NAME', 
        },
        locus_link_id => {
            is => 'Number', 
            column_name => 'LOCUS_LINK_ID', 
        },
    ],
    data_source => 'Genome::DataSource::Oltp',
};

1;

