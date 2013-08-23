package Genome::DruggableGene::DrugAlternateNameReport;

use strict;
use warnings;

use Genome;

class Genome::DruggableGene::DrugAlternateNameReport {
    is => 'UR::Object',
    id_generator => '-uuid',
    table_name => 'dgidb.drug_name_report_association',
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
            constraint_name => 'drug_name_report_association_drug_name_report_id_fkey',
        },
        alternate_name => {is => 'Text'},
        nomenclature => { is => 'Text'},
        description => {
            is => 'Text',
            is_optional => 1,
        },
    ],
    doc => 'Claim regarding an alternate name for a drug name',
};

1;
