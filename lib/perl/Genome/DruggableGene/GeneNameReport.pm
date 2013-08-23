package Genome::DruggableGene::GeneNameReport;

use strict;
use warnings;

use Genome;
use List::MoreUtils qw/ uniq /;

class Genome::DruggableGene::GeneNameReport {
    is => 'UR::Object',
    id_generator => '-uuid',
    table_name => 'dgidb.gene_name_report',
    schema_name => 'dgidb',
    data_source => 'Genome::DataSource::Dgidb',
    id_by => [
        id => {is => 'Text'},
    ],
    has => [
        name => { is => 'Text'},
        nomenclature => { is => 'Text'},
        description => {
            is => 'Text',
            is_optional => 1,
        },
        gene_alt_names => {
            is => 'Genome::DruggableGene::GeneAlternateNameReport',
            reverse_as => 'gene',
            is_many => 1,
        },
        alternate_names => {
            via => 'gene_alt_names',
            to => 'alternate_name',
            is_many => 1,
        },
        gene_categories => {
            is => 'Genome::DruggableGene::GeneCategoryReport',
            reverse_as => 'gene',
            is_many => 1,
        },
        interactions => {
            is => 'Genome::DruggableGene::DrugGeneInteractionReport',
            reverse_as => 'gene',
            is_many => 1,
        },
        drugs => {
            is => 'Genome::DruggableGene::DrugNameReport',
            via => 'interactions',
            to => 'drug',
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
        is_kinase => {
            calculate => q{
                return 1 if grep($_->alternate_name =~ /kinase/i, $self->gene_alt_names);
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

sub original_data_source_url {
    my $self = shift;
    my $base_url = $self->citation->base_url;
    my $source_id = $self->name;
    my $url;
    if($self->source_db_name eq 'DrugBank'){
        $url = join('/', $base_url, 'molecules', $source_id . '?as=target');
    }elsif($self->source_db_name eq 'TTD'){
        $url = $base_url . 'Detail.asp?ID=' . $source_id;
    }elsif($self->source_db_name eq 'Ensembl'){
        $url = $base_url . $source_id;
    }elsif($self->source_db_name eq 'Entrez'){
        $url = $base_url . $source_id;
    }elsif($self->source_db_name eq 'GO'){
        my ($go_short_name_and_id) = map($_->category_value, grep($_->category_name eq 'go_short_name_and_id', $self->gene_categories));
        my (undef, $go_id) = split('_', $go_short_name_and_id);
        $go_id =~ s/^go//;
        $url = $base_url . $go_id;
    }else{
        $url = join('', $base_url, $source_id);
    }

    return $url;
}

sub convert_to_entrez  {
    my $class = shift;
    my @gene_identifiers = @_;
    my ($entrez_gene_symbol_matches, $entrez_id_matches, $ensembl_id_matches, $uniprot_id_matches);
    my $intermediate_genes;

    @gene_identifiers = $class->_strip_version_numbers(@gene_identifiers);

    ($entrez_gene_symbol_matches, @gene_identifiers) = $class->_match_as_entrez_gene_symbol(@gene_identifiers);

    if(@gene_identifiers){
        ($entrez_id_matches, @gene_identifiers) = $class->_match_as_entrez_id(@gene_identifiers);
    }

    if(@gene_identifiers){
        ($ensembl_id_matches, @gene_identifiers) = $class->_match_as_ensembl_id(@gene_identifiers);
    }

    if(@gene_identifiers){
        ($uniprot_id_matches, $intermediate_genes, @gene_identifiers) = $class->_match_as_uniprot_id(@gene_identifiers);
    }

    my $merged_conversion_results = $class->_merge_conversion_results($entrez_gene_symbol_matches, $entrez_id_matches, $ensembl_id_matches, $uniprot_id_matches);

    return $merged_conversion_results, $intermediate_genes, @gene_identifiers;
}

sub _match_as_entrez_gene_symbol {
    my $class = shift;
    my @gene_identifiers = @_;
    my %matched_identifiers;
    my @unmatched_identifiers;

    my @entrez_gene_alt_names = Genome::DruggableGene::GeneAlternateNameReport->get(nomenclature => ['entrez_gene_symbol', 'entrez_gene_synonym'], alternate_name => \@gene_identifiers);
    return {}, @gene_identifiers unless @entrez_gene_alt_names;
    for my $gene_identifier(@gene_identifiers){
        my @associations_for_identifier = grep($_->alternate_name eq $gene_identifier, @entrez_gene_alt_names);
        if(@associations_for_identifier){
            my @genes_for_identifier = map($_->gene, @associations_for_identifier);
            @genes_for_identifier = uniq @genes_for_identifier;
            $matched_identifiers{$gene_identifier} = \@genes_for_identifier;
        }else{
            push @unmatched_identifiers, $gene_identifier;
        }
    }

    return \%matched_identifiers, @unmatched_identifiers;
}

sub _match_as_entrez_id {
    my $class = shift;
    my @gene_identifiers = @_;
    my %matched_identifiers;
    my @unmatched_identifiers;

    my @entrez_genes = Genome::DruggableGene::GeneNameReport->get(nomenclature => 'entrez_id', name => \@gene_identifiers);
    return {}, @gene_identifiers unless @entrez_genes;
    for my $gene_identifier (@gene_identifiers){
        my @reports_for_identifier = grep($_->name eq $gene_identifier, @entrez_genes);
        if(@reports_for_identifier){
            $matched_identifiers{$gene_identifier} = \@reports_for_identifier;

        }else{
            push @unmatched_identifiers, $gene_identifier;
        }
    }

    return \%matched_identifiers, @unmatched_identifiers;
}

sub _match_as_ensembl_id {
    my $class = shift;
    my @gene_identifiers = @_;
    my %matched_identifiers;
    my @unmatched_identifiers;

    my @genes = Genome::DruggableGene::GeneNameReport->get(source_db_name => 'Ensembl', name => \@gene_identifiers);
    for my $gene_identifier(@gene_identifiers){
        my @reports_for_identifier = grep($_->name eq $gene_identifier, @genes);
        unless(@reports_for_identifier){
            push @unmatched_identifiers, $gene_identifier;
            next;
        }
        my @temporary_identifiers = (map($_->name, @reports_for_identifier), map($_->alternate_name, map($_->gene_alt_names, @reports_for_identifier)));
        my ($matched_temporary_identifiers) = $class->_match_as_entrez_gene_symbol(@temporary_identifiers);
        my @complete_reports_for_identifier = map(@{$matched_temporary_identifiers->{$_}}, keys %$matched_temporary_identifiers);
        if(@complete_reports_for_identifier){
            my @complete_reports_for_identifier = uniq @complete_reports_for_identifier;
            $matched_identifiers{$gene_identifier} = \@complete_reports_for_identifier;
        }else{
            push @unmatched_identifiers, $gene_identifier;
        }
    }
    return \%matched_identifiers, @unmatched_identifiers;
}

sub _match_as_uniprot_id {
    my $class = shift;
    my @gene_identifiers = @_;
    my %matched_identifiers;
    my %intermediate_results_for_identifiers;
    my @unmatched_identifiers;

    my @uniprot_associations = Genome::DruggableGene::GeneAlternateNameReport->get(nomenclature => 'uniprot_id', alternate_name => [@gene_identifiers]);
    for my $gene_identifier(@gene_identifiers){
        my @associations_for_identifier = grep($_->alternate_name => @uniprot_associations);
        unless(@associations_for_identifier){
            push @unmatched_identifiers, $gene_identifier;
            next;
        }
        my @uniprot_reports_for_identifier = map($_->gene, @associations_for_identifier);
        @uniprot_reports_for_identifier = uniq @uniprot_reports_for_identifier;
        $intermediate_results_for_identifiers{$gene_identifier} = \@uniprot_reports_for_identifier;
        my @temporary_identifiers = ( (map{$_->name}@uniprot_reports_for_identifier), map{$_->alternate_name} grep{$_->nomenclature ne 'uniprot_id'} map{$_->gene_alt_names} @uniprot_reports_for_identifier);
        my ($matched_temporary_identifiers) = $class->_match_as_entrez_gene_symbol(@temporary_identifiers);
        my @complete_reports_for_identifier = map(@{$matched_temporary_identifiers->{$_}}, keys %$matched_temporary_identifiers);
        @complete_reports_for_identifier = uniq @complete_reports_for_identifier;
        $matched_identifiers{$gene_identifier} = \@complete_reports_for_identifier;
    }

    return \%matched_identifiers, \%intermediate_results_for_identifiers, @unmatched_identifiers;
}

sub _strip_version_numbers{
    my $class = shift;
    my @gene_identifiers = @_;
    my @updated_gene_identifiers;
    #If the incoming gene identifier has a trailing version number, strip it off before comparison
    for my $gene_identifier (@gene_identifiers){
        if ($gene_identifier =~ /(.*)\.\d+$/){
            $gene_identifier = $1;
        }
        push @updated_gene_identifiers, $gene_identifier;
    }
    return @updated_gene_identifiers;
}

sub _merge_conversion_results{
    my $class = shift;
    my @conversion_results = @_;
    my %merged = ();
    for my $conversion_result (@conversion_results){
        %merged = (%merged, %$conversion_result) if $conversion_result;
    }
    return \%merged;
}

sub human_readable_name {
    my $self = shift;
    #Pick a name that is likely human readable
    my ($name) = map{$_->alternate_name}grep{$_->nomenclature eq 'entrez_gene_symbol'} $self->gene_alt_names;
    $name = $self->name if not $name and $self->source_db_name eq 'GO';
    ($name) = $self->alternate_names unless $name;
    $name = $self->name unless $name;
    return $name;
}

1;
