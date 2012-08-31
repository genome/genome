package Genome::DruggableGene::Command::GeneNameReport::LookupInteractions;

use strict;
use warnings;
use Genome;
use List::MoreUtils qw/ uniq /;
use Set::Scalar;

class Genome::DruggableGene::Command::GeneNameReport::LookupInteractions {
    is => 'Genome::Command::Base',
    has_optional => [
        output_file => {
            is => 'Text',
            is_input => 1,
            is_output=> 1,
            doc => "Output interactions to specified file. Defaults to STDOUT if no file is supplied.",
            default => "STDOUT",
        },
        gene_file => {
            is => 'Path',
            is_input => 1,
            doc => 'Path to a list of gene identifiers',
            shell_args_position => 1,
        },
        gene_identifiers => {
            is => 'Text',
            is_many => 1,
            doc => 'Array of gene identifiers',
        },
        filter => {
            is => 'Text',
            doc => 'Filter results based on the parameters.  See below for how to.',
            shell_args_position => 2,
        },
        headers => {
            is => 'Boolean',
            default => 1,
            doc => 'Do include headers',
        },
    ],
    has_transient_optional => [
        output => {
            is => 'Text',
            is_many => 1,
            doc => 'Save output for caller',
        },
        no_match_genes => {
            is_many => 1,
            is => 'Text',
        },
        no_interaction_genes => {
            is_many => 1,
            is => 'Genome::DruggableGene::GeneNameReport',
        },
        filtered_out_interactions => {
            is_many => 1,
            is => 'Genome::DruggableGene::GeneNameReport',
        },
        interactions => {
            is_many => 1,
            is => 'Genome::DruggableGene::DrugGeneInteractionReport',
        },
        identifier_to_genes => {
            is => 'HASH',
        },
        gene_to_identifiers => {
            is => 'HASH',
        }
    ],
};

sub help_brief { 'Lookup drug-gene interactions by gene identifiers' }

sub help_synopsis {
    "genome druggable-gene gene-name-report lookup-interactions --gene-file ./gene_file.txt --filter 'drug.is_withdrawn=0,drug.is_nutraceutical=0,is_potentiator=0,is_untyped=0,drug.is_antineoplastic=1,gene.is_kinase=0'"

}

sub help_detail {
    return <<EOS
LookupInteractions expects a gene-file consisting of a list of gene identifiers one per line and outputs one interaction per line.

Example Syntax:

genome druggable-gene gene-name-report lookup-interactions --gene-file=gene_names.txt --filter='drug.is_approved=1,drug.is_withdrawn=0,drug.is_nutraceutical=0,interaction_attributes.name=is_known_action,interaction_attributes.value=yes,is_potentiator=0'

genome druggable-gene gene-name-report lookup-interactions --gene-file=gene_names.txt --filter='drug.is_withdrawn=0,drug.is_nutraceutical=0,interaction_attributes.name=is_known_action,interaction_attributes.value=yes,is_potentiator=0'

genome druggable-gene gene-name-report lookup-interactions --gene-file=gene_names.txt --filter='drug.is_withdrawn=0,drug.is_nutraceutical=0,is_potentiator=0,(is_untyped=0 or is_known_action=1)';

genome druggable-gene gene-name-report lookup-interactions --gene-file=gene_names.txt --filter='drug.is_withdrawn=0,drug.is_nutraceutical=0,is_potentiator=0,is_inhibitor=1,(is_untyped=0 or is_known_action=1)';

genome druggable-gene gene-name-report lookup-interactions --gene-file=gene_names.txt --filter='drug.is_withdrawn=0,drug.is_nutraceutical=0,is_potentiator=0,gene.is_kinase=1,(is_untyped=0 or is_known_action=1)';

genome druggable-gene gene-name-report lookup-interactions --gene-file=gene_names.txt --filter='drug.is_withdrawn=0,drug.is_nutraceutical=0,is_potentiator=0,drug.is_antineoplastic=1,(is_untyped=0 or is_known_action=1)';
EOS
}

sub execute {
    my $self = shift;

    my @gene_identifiers;
    @gene_identifiers = $self->_read_gene_file();
    @gene_identifiers = $self->gene_identifiers unless @gene_identifiers;
    unless(@gene_identifiers){
        $self->output($self->error_message('No genes found'));
        return;
    }
    my ($unmatched_genes, @genes) = $self->lookup_gene_identifiers(@gene_identifiers);
    my ($no_interaction_genes, $filtered_out_interactions, $interactions) = $self->get_interactions(@genes);
    my %grouped_interactions = $self->group_interactions_by_drug(@$interactions);
    $self->print_grouped_interactions(%grouped_interactions);

    $self->no_match_genes($unmatched_genes);
    $self->no_interaction_genes([@$no_interaction_genes]);
    $self->filtered_out_interactions([@$filtered_out_interactions]);
    $self->interactions([@$interactions]);

    return 1;
}

sub lookup_gene_identifiers {
    my $self = shift;
    my @gene_identifiers = @_;

    my ($entrez_genes, $entrez_intermediate_genes, @unmatched_gene_identifiers) = Genome::DruggableGene::GeneNameReport->convert_to_entrez(@gene_identifiers);
    my %genes = $self->_find_genes_for_identifiers(@gene_identifiers);
    $self->identifier_to_genes(\%genes);
    $self->_create_gene_to_identifiers();

    my @complete_genes;
    for my $gene_identifier (@gene_identifiers){
        my $entrez_genes_for_identifier = $entrez_genes->{$gene_identifier};
        push @complete_genes, @$entrez_genes_for_identifier if $entrez_genes_for_identifier;
        my $entrez_intermediate_genes_for_identifier = $entrez_intermediate_genes->{$gene_identifier};
        push @complete_genes, @$entrez_intermediate_genes_for_identifier if $entrez_intermediate_genes_for_identifier;
        my $genes_for_identifier = $genes{$gene_identifier};
        push @complete_genes, @$genes_for_identifier if $genes_for_identifier;
    }

    return \@unmatched_gene_identifiers, uniq @complete_genes;
}

sub _find_genes_for_identifiers {
    my $self = shift;
    my @gene_identifiers = @_;
    my %results;

    my @genes = Genome::DruggableGene::GeneNameReport->get($self->_chunk_in_clause_list('Genome::DruggableGene::GeneNameReport', 'name', '', @gene_identifiers));
    my @gene_alternates = Genome::DruggableGene::GeneAlternateNameReport->get($self->_chunk_in_clause_list('Genome::DruggableGene::GeneAlternateNameReport', 'alternate_name', '',  @gene_identifiers));
    my @ids = map($_->gene_id, @gene_alternates);
    @ids = uniq @ids;
    Genome::DruggableGene::GeneNameReport->get($self->_chunk_in_clause_list('Genome::DruggableGene::GeneNameReport', 'id', '', @ids));
    push @ids, map($_->id, @genes);
    Genome::DruggableGene::GeneAlternateNameReport->get($self->_chunk_in_clause_list('Genome::DruggableGene::GeneAlternateNameReport', 'gene_id', '', @ids));
    for my $gene_identifier(@gene_identifiers) {
        my @reports_for_identifier = grep($_->name eq $gene_identifier, @genes);
        my @associations_for_identifier = grep($_->alternate_name eq $gene_identifier, @gene_alternates);
        my @report_ids = map($_->gene_id, @associations_for_identifier);
        @reports_for_identifier = (@reports_for_identifier, Genome::DruggableGene::GeneNameReport->get($self->_chunk_in_clause_list('Genome::DruggableGene::GeneNameReport', 'id', '', @report_ids)));
        @reports_for_identifier = uniq @reports_for_identifier;
        $results{$gene_identifier} = \@reports_for_identifier;
    }
    return %results;
}

sub get_interactions {
    my $self = shift;
    my $genes = Set::Scalar->new(@_);

    my @gene_ids = map($_->id, @$genes);
    @gene_ids = uniq @gene_ids;
    my $unfiltered_interactions = Set::Scalar->new (
        Genome::DruggableGene::DrugGeneInteractionReport->get(
            $self->_chunk_in_clause_list('Genome::DruggableGene::DrugGeneInteractionReport', 'gene_id', '', @gene_ids)
        )
    );

    my $genes_in_unfiltered_interactions = Set::Scalar->new(map{$_->gene}@$unfiltered_interactions);
    my $no_interaction_genes = $genes - $genes_in_unfiltered_interactions;

    my @drug_ids = map($_->drug_id, @$unfiltered_interactions);
    Genome::DruggableGene::DrugNameReport->get(\@drug_ids);
    Genome::DruggableGene::DrugCategoryReport->get($self->_chunk_in_clause_list('Genome::DruggableGene::DrugCategoryReport', 'drug_id', '', @drug_ids));
    Genome::DruggableGene::DrugGeneInteractionReportAttribute->get( $self->_chunk_in_clause_list( 'Genome::DruggableGene::DrugGeneInteractionReportAttribute', 'interaction_id', '', map($_->id, @$unfiltered_interactions)));

    my $interactions = Set::Scalar->new(
        Genome::DruggableGene::DrugGeneInteractionReport->get(
            $self->_chunk_in_clause_list('Genome::DruggableGene::DrugGeneInteractionReport', 'gene_id', $self->filter, @gene_ids)
        )
    );

    my $filtered_out_interactions = $unfiltered_interactions - $interactions;

    return $no_interaction_genes, $filtered_out_interactions, $interactions;
}

sub group_interactions_by_drug{
    my $self = shift;
    my @interactions = @_;
    my %grouped_interactions = ();

    for my $interaction (@interactions){
        my $drug_id = $interaction->drug_id;
        if($grouped_interactions{$drug_id}){
            my @temp = @{$grouped_interactions{$drug_id}};
            push @temp, $interaction;
            $grouped_interactions{$drug_id} = \@temp;
        }
        else{
            $grouped_interactions{$drug_id} = [$interaction];
        }
    }

    return %grouped_interactions;
}

sub print_grouped_interactions{
    my $self = shift;
    my %grouped_interactions = @_;

    my $output_file = $self->output_file;
    my $output_fh;
    if ($self->output_file =~ /STDOUT/i) {
        $output_fh = 'STDOUT';
    }else{
        $output_fh = IO::File->new($self->output_file, 'w');
        unless($output_fh){
            $self->error_message("Could not open file " . $self->output_file . " : $@");
            return;
        }
    }

    my @headers = qw/
    drug
    drug_nomenclature
    drug_primary_name
    drug_alternate_names
    drug_brands
    drug_types
    drug_groups
    drug_categories
    drug_source_db_name
    drug_source_db_version
    gene_identifiers
    gene
    gene_nomenclature
    gene_alternate_names
    gene_source_db_name
    gene_source_db_version
    entrez_gene_name
    entrez_gene_synonyms
    interaction_types
    /;
    if($self->headers){
        $output_fh->print(join("\t", @headers), "\n");
        $self->output([join("\t", @headers)]);
    }

    my @drugs = Genome::DruggableGene::DrugNameReport->get($self->_chunk_in_clause_list('Genome::DruggableGene::DrugNameReport', 'id', '', keys %grouped_interactions));
    for my $drug_id (keys %grouped_interactions){
        for my $interaction (@{$grouped_interactions{$drug_id}}){
            $output_fh->print($self->_build_interaction_line($interaction), "\n");
            $self->output([$self->output , $self->_build_interaction_line($interaction)]);
        }
    }

    unless($self->output_file =~ /STDOUT/i){
        $output_fh->close;
    }

    return 1;
}

sub _build_interaction_line {
    my $self = shift;
    my $interaction = shift;
    my $drug = $interaction->drug;
    my $gene = $interaction->gene;
    my $gene_alternate_names = join(':', map($_->alternate_name, $gene->gene_alt_names));
    my $gene_identifiers = join(':', sort @{$self->gene_to_identifiers->{$gene->id}});
    my ($entrez_gene_name, $entrez_gene_synonyms) = $self->_create_entrez_gene_outputs($gene_identifiers);
    my ($drug_primary_name) = map($_->alternate_name, grep($_->nomenclature =~ /primary/i, $drug->drug_alt_names));
    my $drug_alternate_names = join(':', map($_->alternate_name, grep($_->nomenclature !~ /primary/i && $_->nomenclature ne 'drug_brand', $drug->drug_alt_names)));
    my $drug_brands = join(':', map($_->alternate_name, grep($_->nomenclature eq 'drug_brand', $drug->drug_alt_names)));
    my $drug_types = join(';', map($_->category_value, grep($_->category_name eq 'drug_type', $drug->drug_categories)));
    my $drug_groups = join(';', map($_->category_value, grep($_->category_name eq 'drug_group', $drug->drug_categories)));
    my $drug_categories = join(';', map($_->category_value, grep($_->category_name eq 'drug_category', $drug->drug_categories)));
    my $interaction_types = join(':', $interaction->interaction_types);
    my $interaction_line = join("\t", $drug->name, $drug->nomenclature, $drug_primary_name,
        $drug_alternate_names, $drug_brands, $drug_types, $drug_groups, $drug_categories, $drug->source_db_name, $drug->source_db_version,
        $gene_identifiers, $gene->name, $gene->nomenclature, $gene_alternate_names,
        $gene->source_db_name, $gene->source_db_version, $entrez_gene_name, $entrez_gene_synonyms, $interaction_types);
    return $interaction_line;
}

sub _create_entrez_gene_outputs{
    my $self = shift;
    my @gene_identifiers = split(':', shift);
    my $entrez_gene_output = "";
    my $entrez_gene_synonyms_output = "";
    my $entrez_delimiter = '|';
    my $sub_delimiter = '/';
    for my $gene_identifier (@gene_identifiers){
        my @genes = @{$self->identifier_to_genes->{$gene_identifier}};
        my @entrez_genes = grep($_->nomenclature eq 'entrez_id', @genes);
        if(@entrez_genes){
            for my $entrez_gene (sort {$a->name cmp $b->name} @entrez_genes){
                my ($entrez_gene_symbol) = sort map($_->alternate_name, grep($_->nomenclature eq 'entrez_gene_symbol', $entrez_gene->gene_alt_names));
                my @entrez_gene_synonyms = sort map($_->alternate_name, grep($_->nomenclature eq 'entrez_gene_synonym', $entrez_gene->gene_alt_names));
                $entrez_gene_output = $entrez_gene_output . ($entrez_gene_output ? $entrez_delimiter : '') . $entrez_gene_symbol;
                $entrez_gene_synonyms_output = $entrez_gene_synonyms_output . ($entrez_gene_synonyms_output ? $entrez_delimiter : '') . join($sub_delimiter, @entrez_gene_synonyms);
            }
        }else{
            $entrez_gene_output = $entrez_gene_output . $entrez_delimiter;
            $entrez_gene_synonyms_output = $entrez_gene_synonyms_output . $entrez_delimiter;
        }
    }
    return ($entrez_gene_output, $entrez_gene_synonyms_output);
}

sub _read_gene_file{
    my $self = shift;
    my $gene_file = $self->gene_file || return;
    my @gene_identifiers;

    my $gene_fh = Genome::Sys->open_file_for_reading($gene_file);

    while (my $gene_identifier = <$gene_fh>){
        chomp $gene_identifier;
        push @gene_identifiers, $gene_identifier;
    }

    $gene_fh->close;

    unless(@gene_identifiers){
        $self->error_message('No gene identifiers in gene_file ' . $self->gene_file . ', exiting');
        return;
    }

    return @gene_identifiers;
}

sub _chunk_in_clause_list{
    my $self = shift;
    my $target_class = shift;
    my $property_name = shift;
    my $filter = shift;
    my @values = @_;

    unless(@values){
        my $boolexpr = $target_class->define_boolexpr($property_name => []);
        return $boolexpr;
    }

    my @chunked_values = [@values];

    my $filter_bx;
    if($filter) {
        my %extra;
        ($filter_bx, %extra) = UR::BoolExpr->resolve_for_string($target_class, $filter);

        $self->error_message( sprintf('Unrecognized field(s): %s', join(', ', keys %extra)) )
            and return if %extra;
    }

    my $boolexpr;
    my @filter_params = ($filter_bx? $filter_bx->params_list : ());
    if(@filter_params and $filter_params[0] eq '-or') {
        $boolexpr = $target_class->define_boolexpr(
            '-or' => [
                map { [$property_name => \@values, @$_] } @{$filter_params[1]},
            ]
        );
    } else {
        $boolexpr = $target_class->define_boolexpr($property_name => \@values, @filter_params);
    }

    return $boolexpr;
}

sub _create_gene_to_identifiers {
    my $self = shift;
    my %identifier_to_genes = %{$self->identifier_to_genes};
    my %gene_to_identifiers;
    for my $identifier (keys %identifier_to_genes){
        my @genes = @{$identifier_to_genes{$identifier}};
        for my $gene (@genes){
            if($gene_to_identifiers{$gene->id}){
                push @{$gene_to_identifiers{$gene->id}}, $identifier;
            }else{
                $gene_to_identifiers{$gene->id} = [$identifier]
            }
        }
    }

    $self->gene_to_identifiers(\%gene_to_identifiers);
}

1;
