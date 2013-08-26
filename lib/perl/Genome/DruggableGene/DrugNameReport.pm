package Genome::DruggableGene::DrugNameReport;

use strict;
use warnings;

use Genome;

class Genome::DruggableGene::DrugNameReport {
    is => 'UR::Object',
    id_generator => '-uuid',
    table_name => 'dgidb.drug_name_report',
    schema_name => 'dgidb',
    data_source => 'Genome::DataSource::Dgidb',
    id_by => [
        id => { is => 'Text'},
    ],
    has => [
        name => { is => 'Text'},
        nomenclature => { is => 'Text'},
        description => {
            is => 'Text',
            is_optional => 1,
        },
        drug_alt_names => {
            is => 'Genome::DruggableGene::DrugAlternateNameReport',
            reverse_as => 'drug',
            is_many => 1,
        },
        alternate_names => {
            via => 'drug_alt_names',
            to => 'alternate_name',
            is_many => 1,
        },
        drug_categories => {
            is => 'Genome::DruggableGene::DrugCategoryReport',
            reverse_as => 'drug',
            is_many => 1,
        },
        interactions => {
            is => 'Genome::DruggableGene::DrugGeneInteractionReport',
            reverse_as => 'drug',
            is_many => 1,
        },
        genes => {
            is => 'Genome::DruggableGene::GeneNameReport',
            via => 'interactions',
            to => 'gene',
            is_many => 1,
        },
        source_db_name => {
            via => 'citation',
            to => 'source_db_name',
        },
        source_db_version => {
            via => 'citation',
            to => 'source_db_version',
        },
        source_db_url => {
            is => 'Text',
            calculate_from => ['source_db_name'],
            calculate => q| Genome::DruggableGene::Citation->source_db_name_to_url($source_db_name) |,
        },
        citation => {
            is => 'Genome::DruggableGene::Citation',
            id_by => 'citation_id',
        },
        citation_id => {
            is => 'Text',
            implied_by => 'citation',
        },
        is_withdrawn => {
            calculate => q{
                return 1 if grep($_->category_value eq 'withdrawn', $self->drug_categories);
                return 0;
            },
        },
        is_nutraceutical => {
            calculate => q{
                return 1 if grep($_->category_value eq 'nutraceutical', $self->drug_categories);
                return 0;
            },
        },
        is_approved => {
            calculate => q{
                return 1 if grep($_->category_value eq 'approved', $self->drug_categories);
                return 0;
            },
        },
        is_antineoplastic => {
            calculate => q{
                return 1 if grep($_->category_value =~ /antineoplastic/, $self->drug_categories);
                return 0;
            },
        },
    ],
    doc => 'Claim regarding the name of a drug',
};

sub __display_name__ {
    my $self = shift;
    return $self->name . '(' . $self->source_db_name . ' ' . $self->source_db_version . ')';
}

sub source_id {
    my $self = shift;
    my $source_id = $self->name;
    return $source_id;
}

sub original_data_source_url {
    my $self = shift;
    my $base_url = $self->citation->base_url;
    my $source_id = $self->source_id;
    my $url;
    if($self->source_db_name eq 'DrugBank'){
        $url = join('/', $base_url, 'drugs', $source_id);
    }elsif($self->source_db_name eq 'TTD'){
        $url = $base_url . 'DRUG.asp?ID=' . $source_id;
    }else{
        $url = join('', $base_url, $source_id);
    }

    return $url;
}

sub human_readable_name {
    my $self = shift;
    #Pick a name that is likely human readable
    my ($name) =  map{$_->alternate_name}grep{$_->nomenclature eq 'TTD_primary_drug_name' or $_->nomenclature eq 'Primary DrugBank drug name'} $self->drug_alt_names;
    ($name) =  map{$_->alternate_name} $self->drug_alt_names unless $name;
    $name = $self->name unless $name;
    return $name;
}

1;
