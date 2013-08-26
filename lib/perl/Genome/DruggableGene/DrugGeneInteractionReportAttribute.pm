package Genome::DruggableGene::DrugGeneInteractionReportAttribute;

use strict;
use warnings;

use Genome;

class Genome::DruggableGene::DrugGeneInteractionReportAttribute {
    is => 'UR::Object',
    id_generator => '-uuid',
    table_name => 'dgidb.drug_gene_interaction_report_attribute',
    id_by => [
        id => { is => 'Text' },
    ],
    has => [
        interaction_id               => { is => 'Text' },
        drug_gene_interaction_report => { is => 'Genome::DruggableGene::DrugGeneInteractionReport', id_by => 'interaction_id', constraint_name => 'drug_gene_interaction_report_attribute_interaction_id_fkey' },
        name                         => { is => 'Text' },
        value                        => { is => 'Text' },
    ],
    schema_name => 'dgidb',
    data_source => 'Genome::DataSource::Dgidb',
    doc => 'Claim regarding an attribute of a drug gene interaction claim',
};

1;
