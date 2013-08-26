package Genome::DruggableGene::GeneCategoryReport;

use strict;
use warnings;

use Genome;

class Genome::DruggableGene::GeneCategoryReport {
    is => 'UR::Object',
    id_generator => '-uuid',
    table_name => 'dgidb.gene_name_report_category_association',
    schema_name => 'dgidb',
    data_source => 'Genome::DataSource::Dgidb',
    id_by => [
        id => {is => 'Text'},
    ],
    has => [
        gene_id => { is => 'Text', column_name => 'gene_name_report_id'},
        gene => {
            is => 'Genome::DruggableGene::GeneNameReport',
            id_by => 'gene_id',
            constraint_name => 'gene_name_report_category_association_gene_name_report_id_fkey',
        },
        category_name => { is => 'Text' },
        category_value => { is => 'Text' },
        description => {
            is => 'Text',
            is_optional => 1,
        },
    ],
    doc => 'Claim regarding a categorization for a gene name',
};

1;
