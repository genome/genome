package Genome::DruggableGene::GeneNameGroupBridge;

use strict;
use warnings;

use Genome;

class Genome::DruggableGene::GeneNameGroupBridge {
    is => 'UR::Object',
    table_name => 'dgidb.gene_name_group_bridge',
    schema_name => 'dgidb',
    data_source => 'Genome::DataSource::Dgidb',

    id_by => [
        group_id => { is => 'Text', column_name => 'gene_name_group_id'},
        gene_id => { is => 'Text', column_name => 'gene_name_report_id'},
    ],
    has => [
        group => {
            is => 'Genome::DruggableGene::GeneNameGroup',
            id_by => 'group_id',
        },
        gene => {
            is => 'Genome::DruggableGene::GeneNameReport',
            id_by => 'gene_id',
        },
    ],
    doc => 'Associate a gene that is likely synonymous with other genes in this group',
};

1;
