package Genome::DruggableGene::DrugGeneInteractionReport;

use strict;
use warnings;

use Genome;

class Genome::DruggableGene::DrugGeneInteractionReport {
    is => 'UR::Object',
    id_generator => '-uuid',
    table_name => 'dgidb.drug_gene_interaction_report',
    schema_name => 'dgidb',
    data_source => 'Genome::DataSource::Dgidb',
    id_by => [
        id => { is => 'Text' },
    ],
    has => [
        drug_id => { is => 'Text', column_name => 'drug_name_report_id'},
        drug => {
            is => 'Genome::DruggableGene::DrugNameReport',
            id_by => 'drug_id',
            constraint_name => 'drug_gene_interaction_report_drug_name_report_id_fkey',
        },
        drug_name => {
            via => 'drug',
            to => 'name',
        },
        human_readable_drug_name => {
            is => 'text',
            calculate_from => ['drug'],
            calculate => q|
                return $drug->human_readable_name;
            |,
        },
        gene_id => { is => 'Text', column_name => 'gene_name_report_id'},
        gene => {
            is => 'Genome::DruggableGene::GeneNameReport',
            id_by => 'gene_id',
            constraint_name => 'drug_gene_interaction_report_gene_name_report_id_fkey',
        },
        gene_name => {
            via => 'gene',
            to => 'name',
        },
        source_db_name => {
            via => 'citation',
            to => 'source_db_name',
        },
        source_db_version => {
            via => 'citation',
            to => 'source_db_version',
        },
        description => { is => 'Text', is_optional => 1 },
        interaction_attributes => {
            is => 'Genome::DruggableGene::DrugGeneInteractionReportAttribute',
            reverse_as => 'drug_gene_interaction_report',
            is_many => 1,
        },
        citation => {
            is => 'Genome::DruggableGene::Citation',
            id_by => 'citation_id',
        },
        citation_id => {
            is => 'Text',
            implied_by => 'citation',
        },
        gene_group_name => {
            is => 'text',
            calculate_from => ['gene_id'],
            calculate => q|
                my $bridge = Genome::DruggableGene::GeneNameGroupBridge->get(gene_id => $gene_id);
                return $bridge->group->name;
            |,
        },
        is_known_action => {
            calculate => q{
                return 1 if grep($_->name eq 'is_known_action' && $_->value eq 'yes', $self->interaction_attributes);
                return 0;
            },
        },
        interaction_types => {
            via => 'interaction_attributes',
            to => 'value',
            where => [name => 'interaction_type'],
            is_many => 1,
            is_optional => 1,
        },
        is_potentiator => {
            calculate => q|
                my @potentiator = grep($_ =~ /potentiator/, $self->interaction_types);
                return 1 if @potentiator;
                return 0;
            |,
        },
        is_untyped => {
            calculate => q|
                my @na = grep($_ =~ /^na$/, $self->interaction_types);
                return 1 if @na;
                return 0;
            |,
        },
        is_inhibitor => {
            calculate => q|
                my @inhibitor = grep($_ =~ /inhibitor/, $self->interaction_types);
                return 1 if @inhibitor;
                return 0;
            |,
        }
    ],
    doc => 'Claim regarding an interaction between a drug name and a gene name',
};

sub __display_name__ {
    my $self = shift;
    return $self->drug->human_readable_name. ' as ' .  join(' and ',$self->interaction_types) .  ' for ' . $self->gene->human_readable_name;
}

1;
