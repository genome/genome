package Genome::DruggableGene::GeneNameReport::Set::View::Interaction::Xml;

use strict;
use warnings;
use Genome;
use Data::Dumper;
use XML::LibXML;

class Genome::DruggableGene::GeneNameReport::Set::View::Interaction::Xml {
    is => 'Genome::View::Status::Xml',
    has => {
        perspective => { is => 'Text', value => 'interaction' },
    },
    has_optional => [
        data => { is => 'HASH' },
    ],
};

sub _generate_content {
    my $self = shift;

    #create the XML doc and add it to the object
    my $doc = XML::LibXML->createDocument();
    $self->_xml_doc($doc);
    my $drug_gene_interaction = $doc->createElement("drug_gene_interaction");

    my $data = $self->data;

    my ($interactions_node, $filtered_interactions_node) = $self->get_interactions($data->{definite_groups});
    $drug_gene_interaction->addChild($interactions_node);
    $drug_gene_interaction->addChild($filtered_interactions_node);

    ($interactions_node, $filtered_interactions_node) = $self->get_ambiguous_interactions($data->{ambiguous_search_terms});
    $drug_gene_interaction->addChild($interactions_node);
    $drug_gene_interaction->addChild($filtered_interactions_node);
    $drug_gene_interaction->addChild($self->get_missing_interactions($data->{definite_groups}));
    $drug_gene_interaction->addChild($self->get_missing_ambiguous_interactions($data->{ambiguous_search_terms}));
    $drug_gene_interaction->addChild($self->get_search_terms_without_groups($data->{search_terms_without_groups}));

    $doc->setDocumentElement($drug_gene_interaction);
    return $doc->toString(1);
}

sub get_interactions {
    my $self = shift;
    my $groups = shift;
    my $doc = $self->_xml_doc;
    my $interactions_node = $doc->createElement("interactions");
    my $filtered_interactions_node = $doc->createElement("filtered_interactions");

    while (my ($group_name, $group_data) = each %{$groups}){
        my $group = $group_data->{group};
        if($group_data->{interactions}){
            for my $interaction (@{$group_data->{interactions}}){
                $interactions_node->addChild($self->build_interaction_node(
                        $interaction,
                        $group_name,
                        $group_data->{search_terms},
                    ));
            }
        }
        if($group_data->{filtered_interactions}){
            for my $interaction (@{$group_data->{filtered_interactions}}){
                $filtered_interactions_node->addChild($self->build_interaction_node(
                        $interaction,
                        $group_name,
                        $group_data->{search_terms},
                    ));
            }
        }
    }
    return $interactions_node, $filtered_interactions_node;
}

sub get_ambiguous_interactions {
    my $self = shift;
    my $ambiguous_terms_to_gene_groups = shift;
    my $doc = $self->_xml_doc;
    my $interactions_node = $doc->createElement("ambiguous_interactions");
    my $filtered_interactions_node = $doc->createElement("filtered_ambiguous_interactions");

    while (my ($ambiguous_term, $gene_groups) = each %{$ambiguous_terms_to_gene_groups}){
        while (my ($gene_group_name, $gene_group_data) = each %{$gene_groups}){
            my $group = $gene_group_data->{group};

            if($gene_group_data->{interactions}){
                for my $interaction (@{$gene_group_data->{interactions}}){
                    my $interaction_node = ($self->build_interaction_node(
                            $interaction,
                            $gene_group_name,
                            [$ambiguous_term],
                        ));

                    my $matches_node = $doc->createElement('number_of_matches');
                    $matches_node->addChild($doc->createTextNode($gene_group_data->{number_of_matches}));
                    $interaction_node->addChild($matches_node);
                    $interactions_node->addChild($interaction_node);
                }
            }
            if($gene_group_data->{filtered_interactions}){
                for my $interaction (@{$gene_group_data->{filtered_interactions}}){
                    my $interaction_node = ($self->build_interaction_node(
                            $interaction,
                            $gene_group_name,
                            [$ambiguous_term],
                        ));

                    my $matches_node = $doc->createElement('number_of_matches');
                    $matches_node->addChild($doc->createTextNode($gene_group_data->{number_of_matches}));
                    $interaction_node->addChild($matches_node);
                    $filtered_interactions_node->addChild($interaction_node);
                }
            }
        }
    }
    return $interactions_node, $filtered_interactions_node;
}

sub get_missing_interactions {
    my $self = shift;
    my $groups = shift;
    my $doc = $self->_xml_doc;
    my $missing_interactions_node = $doc->createElement("missing_interactions");

    while (my ($group_name, $group_data) = each %{$groups}){
        my $group = $group_data->{group};
        my @search_terms = @{$group_data->{search_terms}};

        if (0 == scalar map{$_->interactions}$group->genes){
            my $item = $doc->createElement('item');
            $missing_interactions_node->addChild($item);
            my $group_node = $doc->createElement('group');
            $group_node->addChild($doc->createTextNode($group_name));
            $item->addChild($group_node);
            my $search_terms_node = $doc->createElement('search_terms');
            $search_terms_node->addChild($doc->createTextNode(join(', ', @search_terms)));
            $item->addChild($search_terms_node);
        }
    }
    return $missing_interactions_node;
}

sub get_missing_ambiguous_interactions {
    my $self = shift;
    my $ambiguous_terms_to_gene_groups = shift;
    my $doc = $self->_xml_doc;
    my $missing_interactions_node = $doc->createElement("missing_ambiguous_interactions");

    while (my ($ambiguous_term, $gene_groups) = each %{$ambiguous_terms_to_gene_groups}){
        while (my ($gene_group_name, $gene_group_data) = each %{$gene_groups}){
            my $group = $gene_group_data->{group};

            if (0 == scalar map{$_->interactions}$group->genes){
                my $item = $doc->createElement('item');
                $missing_interactions_node->addChild($item);
                my $group_node = $doc->createElement('group');
                $group_node->addChild($doc->createTextNode($gene_group_name));
                $item->addChild($group_node);
                my $search_terms_node = $doc->createElement('search_terms');
                $search_terms_node->addChild($doc->createTextNode($ambiguous_term));
                $item->addChild($search_terms_node);
                my $matches_node = $doc->createElement('number_of_matches');
                $matches_node->addChild($doc->createTextNode($gene_group_data->{number_of_matches}));
                $item->addChild($matches_node);
            }
        }
    }
    return $missing_interactions_node;
}

sub get_search_terms_without_groups {
    my $self = shift;
    my $search_terms_without_groups = shift;
    my $doc = $self->_xml_doc;
    my $search_terms_without_groups_node = $doc->createElement("search_terms_without_groups");

    for my $name (@{$search_terms_without_groups}){
        my $item = $doc->createElement('item');
        $item->addChild($doc->createTextNode($name));
        $search_terms_without_groups_node->addChild($item);
    }

    return $search_terms_without_groups_node;
}

sub build_interaction_node {
    my $self = shift;
    my $interaction = shift;
    my $gene_group_name = shift;
    my $search_terms = shift;
    my $doc = $self->_xml_doc;
    my $drug_name = $interaction->drug_name,
    my $human_readable_drug_name = $interaction->drug->human_readable_name,
    my $gene_name = $interaction->gene_name,
    my $interaction_types = [$interaction->interaction_types],

    my $item = $doc->createElement('item');
    my $drug_node = $doc->createElement('drug');
    $drug_node->addChild($doc->createAttribute('key', 'drug_name'));
    $drug_node->addChild($doc->createTextNode($drug_name));
    $item->addChild($drug_node);
    my $human_readable_drug_name_node = $doc->createElement('human_readable_drug_name');#Eventually replace with drug groups
    $human_readable_drug_name_node->addChild($doc->createTextNode($human_readable_drug_name));
    $item->addChild($human_readable_drug_name_node);
    my $gene_node = $doc->createElement('gene');
    $gene_node->addChild($doc->createAttribute('key', 'gene_name'));
    $gene_node->addChild($doc->createTextNode($gene_name));
    $item->addChild($gene_node);
    my $group_node = $doc->createElement('group');
    $group_node->addChild($doc->createTextNode($gene_group_name));
    $item->addChild($group_node);
    my $interaction_types_node = $doc->createElement('interaction_type');
    $interaction_types_node->addChild($doc->createTextNode(join(', ', @{$interaction_types})));
    $item->addChild($interaction_types_node);
    my $search_terms_node = $doc->createElement('search_terms');
    $search_terms_node->addChild($doc->createTextNode(join(', ', @{$search_terms})));
    $item->addChild($search_terms_node);
    my $source_node = $doc->createElement('source');
    $source_node->addChild($doc->createTextNode($interaction->source_db_name));
    $item->addChild($source_node);
    my $source_url_node = $doc->createElement('source_url');
    $source_url_node->addChild($doc->createTextNode(Genome::DruggableGene::Citation->source_db_name_to_url($interaction->source_db_name)));
    $item->addChild($source_url_node);

    return $item;
}
1;
