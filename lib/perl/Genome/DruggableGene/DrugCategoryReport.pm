package Genome::DruggableGene::DrugCategoryReport;

use strict;
use warnings;

use Genome;

class Genome::DruggableGene::DrugCategoryReport {
    is => 'UR::Object',
    id_generator => '-uuid',
    table_name => 'dgidb.drug_name_report_category_association',
    schema_name => 'dgidb',
    data_source => 'Genome::DataSource::Dgidb',
    id_by => [
        id => {is => 'Text'},
    ],
    has => [
        drug_id => { is => 'Text', column_name => 'drug_name_report_id'},
        drug => {
            is => 'Genome::DruggableGene::DrugNameReport',
            id_by => 'drug_id',
            constraint_name => 'drug_name_report_category_association_drug_name_report_id_fkey',
        },
        category_name => { is => 'Text' },
        category_value => { is => 'Text' },
        description => {
            is => 'Text',
            is_optional => 1,
        },
    ],
    doc => 'Claim regarding categorization of a drug name',
};

1;
