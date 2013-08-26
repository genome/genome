package Genome::DruggableGene::GeneAlternateNameReport;

use strict;
use warnings;

use Genome;

class Genome::DruggableGene::GeneAlternateNameReport {
    is => 'UR::Object',
    id_generator => '-uuid',
    table_name => 'dgidb.gene_name_report_association',
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
            constraint_name => 'gene_name_report_association_gene_name_report_id_fkey',
        },
        alternate_name => {is => 'Text'},
        nomenclature => { is => 'Text'},
        description => {
            is => 'Text',
            is_optional => 1,
        },
    ],
    doc => 'Claim regarding an alternate name for a gene name',
};

1;
