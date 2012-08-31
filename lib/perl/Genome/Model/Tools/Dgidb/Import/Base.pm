package Genome::Model::Tools::Dgidb::Import::Base;

use strict;
use warnings;

use Genome;

class Genome::Model::Tools::Dgidb::Import::Base {
    is => 'Command::V2',
    is_abstract => 1,
    has => [
        version => {
            is => 'Text',
            is_input => 1,
            doc => 'Version identifier for the infile (ex 3)',
        },
        citation_base_url => {
            is => 'Text',
            is_input => 1,
            doc => 'base url string for the citation object',
        },
        citation_site_url => {
            is => 'Text',
            is_input => 1,
            doc => 'site url for the citation object',
        },
        citation_text => {
            is => 'Text',
            is_input => 1,
            doc => 'citation text for this datasource',
        },
    ],
    doc => 'Base class for importing datasets into DGI:DB',
};

sub _create_citation {
    my $self = shift;
    my $source_db_name = shift;
    my $source_db_version = shift;
    my $base_url = shift;
    my $site_url = shift;
    my $citation_text = shift;
    return Genome::DruggableGene::Citation->create(
        source_db_name => $source_db_name,
        source_db_version => $source_db_version,
        base_url => $base_url,
        site_url => $site_url,
        citation => $citation_text,
    );
}

sub _create_drug_name_report {
    my $self = shift;
    my ($name, $citation, $nomenclature, $description) = @_;
    my %params = (
        name => uc $name,
        nomenclature => $nomenclature,
        citation => $citation,
        description => $description,
    );

    my $drug_name_report = Genome::DruggableGene::DrugNameReport->get(%params);
    return $drug_name_report if $drug_name_report;
    return Genome::DruggableGene::DrugNameReport->create(%params);
}

sub _create_drug_alternate_name_report {
    my $self = shift;
    my ($drug_name_report, $alternate_name, $nomenclature, $description) = @_;
    my %params = (
        drug_id => $drug_name_report->id,
        alternate_name => uc $alternate_name,
        nomenclature => $nomenclature,
        description => $description,
    );

    my $drug_alternate_name_report = Genome::DruggableGene::DrugAlternateNameReport->get(%params);
    return $drug_alternate_name_report if $drug_alternate_name_report;
    return Genome::DruggableGene::DrugAlternateNameReport->create(%params);
}

sub _create_drug_category_report {
    my $self = shift;
    my ($drug_name_report, $category_name, $category_value, $description) = @_;
    my %params = (
        drug_id => $drug_name_report->id,
        category_name => $category_name,
        category_value => lc $category_value,
        description => $description,
    );
    my $drug_category_report = Genome::DruggableGene::DrugCategoryReport->get(%params);
    return $drug_category_report if $drug_category_report;
    return Genome::DruggableGene::DrugCategoryReport->create(%params);
}

sub _create_gene_name_report {
    my $self = shift;
    my ($name, $citation, $nomenclature, $description) = @_;
    my %params = (
        name => uc $name,
        nomenclature => $nomenclature,
        citation => $citation,
        description => $description,
    );

    if($name ne 'NA'){
        my $gene_name_report = Genome::DruggableGene::GeneNameReport->get(%params);
        return $gene_name_report if $gene_name_report;
    }
    return Genome::DruggableGene::GeneNameReport->create(%params);
}

sub _create_gene_alternate_name_report {
    my $self = shift;
    my ($gene_name_report, $alternate_name, $nomenclature, $description) = @_;
    my %params = (
        gene_id => $gene_name_report->id,
        alternate_name => uc $alternate_name,
        nomenclature => $nomenclature,
        description => $description,
    );
    my $gene_alternate_name_report = Genome::DruggableGene::GeneAlternateNameReport->get(%params);
    return $gene_alternate_name_report if $gene_alternate_name_report;
    return Genome::DruggableGene::GeneAlternateNameReport->create(%params);
}

sub _create_gene_category_report {
    my $self = shift;
    my ($gene_name_report, $category_name, $category_value, $description) = @_;
    my %params = (
        gene_id => $gene_name_report->id,
        category_name => $category_name,
        category_value => lc $category_value,
        description => $description,
    );
    my $gene_category_report = Genome::DruggableGene::GeneCategoryReport->get(%params);
    return $gene_category_report if $gene_category_report;
    return Genome::DruggableGene::GeneCategoryReport->create(%params);
}

sub _create_interaction_report {
    my $self = shift;
    my ($citation, $drug_name_report, $gene_name_report, $description) = @_;
    my %params = (
        gene_id => $gene_name_report->id,
        drug_id => $drug_name_report->id,
        citation => $citation,
        description =>  $description,
    );

    my $interaction = Genome::DruggableGene::DrugGeneInteractionReport->get(%params);
    return $interaction if $interaction;
    return Genome::DruggableGene::DrugGeneInteractionReport->create(%params);
}

sub _create_interaction_report_attribute {
    my $self = shift;
    my ($interaction, $name, $value) = @_;
    my %params = (
        drug_gene_interaction_report => $interaction,
        name => $name,
        value => lc $value,
    );
    my $attribute = Genome::DruggableGene::DrugGeneInteractionReportAttribute->get(%params);
    return $attribute if $attribute;
    return Genome::DruggableGene::DrugGeneInteractionReportAttribute->create(%params);
}

sub _destroy_and_rebuild_pubchem_and_drug_groups {
    my $citation_text = Genome::DruggableGene::Citation->get(source_db_name => 'PubChem')->citation;
    #nuke drug groups
    Genome::DruggableGene::Command::DrugNameGroup::RemoveAll->execute();
    #nuke pubchem
    Genome::DruggableGene::Command::RemoveDatasource->execute(source_db_name => 'PubChem');
    #reimport pubchem
    Genome::Model::Tools::Dgidb::Import::Pubchem->execute(postgres_host => 'postgres', scratch_dir => '/tmp/', citation_text => $citation_text, generate_uuids_locally => 1); #TODO: handle different postgres hosts
    #regroup_drugs
    Genome::DruggableGene::Command::DrugNameGroup::Generate->execute();
    return 1;
}

1;
