package Genome::DruggableGene::DrugNameGroupBridge;

use strict;
use warnings;

use Genome;

class Genome::DruggableGene::DrugNameGroupBridge {
    is => 'UR::Object',
    table_name => 'dgidb.drug_name_group_bridge',
    schema_name => 'dgidb',
    data_source => 'Genome::DataSource::Dgidb',

    id_by => [
        group_id => { is => 'Text', column_name => 'drug_name_group_id'},
        drug_id => { is => 'Text', column_name => 'drug_name_report_id'},
    ],
    has => [
        group => {
            is => 'Genome::DruggableGene::DrugNameGroup',
            id_by => 'group_id',
        },
        drug => {
            is => 'Genome::DruggableGene::DrugNameReport',
            id_by => 'drug_id',
        },
    ],
    doc => 'Associate a drug that is likely synonymous with other drugs in this group',
};

1;
