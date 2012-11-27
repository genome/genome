package Genome::DruggableGene::Command::RemoveDatasource;

use strict;
use warnings;
use List::MoreUtils qw/ uniq /;
use Genome;

class Genome::DruggableGene::Command::RemoveDatasource {
    is => 'Genome::Command::Base',
    has => [
        source_db_name => {
            is => 'Text',
            doc => 'source_db_name of the datasource to remove (e.g., dGene, DrugBank, Ensembl, Entrez, GO, HopkinsGroom, PharmGKB, RussLampel, TALC, TTD)',
        },
        source_db_version => {
            is => 'Text',
            is_optional => 1,
            doc => 'source_db_version of the datasource to remove.  If unspecified, all versions of the source will be deleted',
        },
    ],
};

sub help_brief { 'Completely remove a DGIDB data source' }

sub help_synopsis {
    return <<HELP
genome druggable-gene remove-datasource --source-db-name=Entrez
HELP
}

sub help_detail { help_brief() }

sub execute {
    my $self = shift;

    my @citation = $self->_fetch_citation;
    my @drugs = $self->_fetch_drugs;
    my @genes = $self->_fetch_genes;
    my @interactions;
    
    for my $gene (@genes){
        for my $a ($gene->gene_alt_names){
            $a->delete;
        }
        for my $c ($gene->gene_categories){
            $c->delete;
        }
        for my $b (Genome::DruggableGene::GeneNameGroupBridge->get(gene_id => $gene->id)){
            $b->delete;
        }
        push @interactions, $gene->interactions;
        $gene->delete;
    }

    for my $drug (@drugs){
        for my $b ($drug->drug_alt_names){
            $b->delete;
        }
        for my $g ($drug->drug_categories){
            $g->delete;
        }
        for my $b (Genome::DruggableGene::DrugNameGroupBridge->get(drug_id => $drug->id)){
            $b->delete;
        }
        push @interactions, $drug->interactions;
        $drug->delete;
    }

    @interactions = uniq @interactions;

    for my $interaction (@interactions){
        for my $att ($interaction->interaction_attributes){
            $att->delete;
        }
        $interaction->delete;
    }

    for my $citation (@citation){
        $citation->delete;
    }

    return 1;
}

sub _fetch_citation {
    my $self = shift;
    my $source_db_name = $self->source_db_name;
    my $source_db_version = $self->source_db_version;
    my @citation; 

    if($source_db_version){
        @citation = Genome::DruggableGene::Citation->get(source_db_name => $source_db_name, source_db_version => $source_db_version);
    }else{
        @citation = Genome::DruggableGene::Citation->get(source_db_name => $source_db_name);
    }
    
    return @citation;
}

sub _fetch_drugs {
    my $self = shift;
    my $source_db_name = $self->source_db_name;
    my $source_db_version = $self->source_db_version;
    my @drugs;

    if($source_db_version){
        @drugs = Genome::DruggableGene::DrugNameReport->get(source_db_name => $source_db_name, source_db_version => $source_db_version);
    }else{
        @drugs = Genome::DruggableGene::DrugNameReport->get(source_db_name => $source_db_name);
    }

    return @drugs;
}

sub _fetch_genes {
    my $self = shift;
    my $source_db_name = $self->source_db_name;
    my $source_db_version = $self->source_db_version;
    my @genes;

    if($source_db_version){
        @genes = Genome::DruggableGene::GeneNameReport->get(source_db_name => $source_db_name, source_db_version => $source_db_version);
    }else{
        @genes = Genome::DruggableGene::GeneNameReport->get(source_db_name => $source_db_name);
    }

    return @genes;
}

1;
