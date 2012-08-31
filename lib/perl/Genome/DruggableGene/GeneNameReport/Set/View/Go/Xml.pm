package Genome::DruggableGene::GeneNameReport::Set::View::Go::Xml;

use strict;
use warnings;
use Genome;
use Data::Dumper;
use XML::LibXML;
use List::MoreUtils qw/ uniq /;

class Genome::DruggableGene::GeneNameReport::Set::View::Go::Xml {
    is => 'Genome::View::Status::Xml',
    has => {
        perspective => { is => 'Text', value => 'go' },
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
    my $go_results = $doc->createElement("go_results");

    my $data = $self->data;
    $go_results->addChild($self->get_go_results('definite_go_results', $data->{definite_groups}));
    $go_results->addChild($self->get_go_results('filtered_definite_go_results', $data->{filtered_definite_groups}));
    $go_results->addChild($self->get_go_summary($data->{definite_groups}));

    $doc->setDocumentElement($go_results);
    return $doc->toString(1);
}

sub get_go_results {
    my $self = shift;
    my $node_name = shift;
    my $groups = shift;
    my $doc = $self->_xml_doc;
    my $go_results_node = $doc->createElement($node_name);

    my %existing_entries;

    while (my ($group_name, $group_data) = each %{$groups}){
        my $group = $group_data->{group};
        my @search_terms = @{$group_data->{search_terms}};
        my @go_genes = grep($_->nomenclature eq 'go_gene_name', $group->genes);
        my @human_readable_names = map($_->alternate_name, grep($_->nomenclature eq 'human_readable_name', map($_->gene_alt_names, @go_genes)));
        for my $human_readable_name (@human_readable_names){
            #skip duplicate entries
            my $entry_key = join(":", $human_readable_name, $group_name, join(', ', @{$group_data->{search_terms}}));
            unless($existing_entries{$entry_key}){
                $go_results_node->addChild($self->build_go_results_node(
                    $human_readable_name,
                    $group_name,
                    $group_data->{search_terms},
                ));
                $existing_entries{$entry_key}++;
            }
        }
    }

    return $go_results_node;
}

sub get_go_summary {
    my $self = shift;
    my $groups = shift;
    my $doc = $self->_xml_doc;
    my $go_summary_node = $doc->createElement("go_results_summary");
    my %go_summary;

    while (my ($group_name, $group_data) = each %{$groups}){
        my $group = $group_data->{group};
        my @go_genes = grep($_->nomenclature eq 'go_gene_name', $group->genes);
        my @human_readable_names = map($_->alternate_name, grep($_->nomenclature eq 'human_readable_name', map($_->gene_alt_names, @go_genes)));
        for my $human_readable_name (@human_readable_names){
            $go_summary{$human_readable_name} = ($go_summary{$human_readable_name} ? join(",", $go_summary{$human_readable_name}, $group_name) : $group_name);
        }
    }

    while(my ($human_readable_name, $group_names) = each %go_summary){
        my @group_names = split(',', $group_names);
        my @uniq_group_names = uniq @group_names;
        $go_summary_node->addChild($self->build_go_summary_node(
            $human_readable_name,
            scalar(@uniq_group_names),
            join(',', @uniq_group_names)
        ));
    }

    return $go_summary_node;
}

sub build_go_results_node {
    my $self = shift;
    my $go_category_name = shift;
    my $gene_group_name = shift;
    my $search_terms = shift;
    my $doc = $self->_xml_doc;

    my $item = $doc->createElement('item');

    my $gene_group_name_node = $doc->createElement('gene_group_name');
    $gene_group_name_node->addChild($doc->createTextNode($gene_group_name));
    $item->addChild($gene_group_name_node);
    my $category_name_node = $doc->createElement('category_name');
    $category_name_node->addChild($doc->createAttribute('key', 'category_name'));
    $category_name_node->addChild($doc->createTextNode($go_category_name));
    $item->addChild($category_name_node);
    my $search_terms_node = $doc->createElement('search_terms');
    $search_terms_node->addChild($doc->createTextNode(join(', ', @{$search_terms})));
    $item->addChild($search_terms_node);

    return $item;
}

sub build_go_summary_node {
    my $self = shift;
    my $category_name = shift;
    my $count = shift;
    my $group_names = shift;
    my $doc = $self->_xml_doc;
    my $item = $doc->createElement('item');

    my $go_category_name_node = $doc->createElement('go_category_name');
    $go_category_name_node->addChild($doc->createTextNode($category_name));
    $item->addChild($go_category_name_node);
    my $count_node = $doc->createElement('gene_count');
    $count_node->addChild($doc->createTextNode($count));
    $item->addChild($count_node);
    my $gene_names_node = $doc->createElement('gene_names');
    $gene_names_node->addChild($doc->createTextNode($group_names));
    $item->addChild($gene_names_node);

    return $item;
}

1;
