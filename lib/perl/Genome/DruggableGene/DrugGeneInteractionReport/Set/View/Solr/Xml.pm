package Genome::DruggableGene::DrugGeneInteractionReport::Set::View::Solr::Xml;

use strict;
use warnings;

use Genome;

class Genome::DruggableGene::DrugGeneInteractionReport::Set::View::Solr::Xml {
    is => 'Genome::View::Solr::Xml',
    has_field => [
        type => {
            is => 'Text',
            default => 'drug-gene-interaction'
        },
        display_type => {
            is  => 'Text',
            default => 'DrugGeneInteraction',
        },
        display_icon_url => {
            is  => 'Text',
            default => 'genome_druggable-gene_drug-gene-interaction_32',
        },
        display_url0 => {
            is => 'Text',
            calculate_from => ['subject'],
            calculate => q{
                return '/view/genome/druggable-gene/drug-gene-interaction-report/set/status.html?drug_name=' . ($subject->members)[0]->drug_name() . '&gene_name=' . ($subject->members)[0]->gene_name();
            },
        },
        display_label1 => {
            is  => 'Text',
        },
        display_url1 => {
            is  => 'Text',
        },
        display_label2 => {
            is  => 'Text',
        },
        display_url2 => {
            is  => 'Text',
        },
        display_label3 => {
            is  => 'Text',
        },
        display_url3 => {
            is  => 'Text',
        },
        display_title => {
            is => 'Text',
            calculate_from => ['subject'],
            calculate => q{ ($subject->members)[0]->__display_name__ },
        },
        title => {
            is => 'Text',
            calculate_from => ['subject'],
            calculate => q{
                join(' ',($subject->members)[0]->interaction_types) .
                ' ' .  ($subject->members)[0]->drug_name .
                ' ' .  ($subject->members)[0]->gene_name .
                ' druggablegene'
            },
        },
        default_aspects => {
            is => 'ARRAY',
            default => [
                {
                    name => 'nomenclature',
                    position => 'content',
                },
                {
                    name => 'source_db_name',
                    position => 'content',
                },
                {
                    name => 'source_db_version',
                    position => 'content',
                },
            ],
        },
    ],
};

1;
