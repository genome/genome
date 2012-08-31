package Genome::DruggableGene::GeneNameReport::Set::View::Solr::Xml;

use strict;
use warnings;

use Genome;

class Genome::DruggableGene::GeneNameReport::Set::View::Solr::Xml {
    is => 'Genome::View::Solr::Xml',
    has_field => [
        type => {
            is => 'Text',
            default => 'gene-name'
        },
        display_type => {
            is  => 'Text',
            default => 'GeneName',
        },
        display_icon_url => {
            is  => 'Text',
            default => 'genome_druggable-gene_gene-name_32',
        },
        display_url0 => {
            is => 'Text',
            calculate_from => ['subject'],
            calculate => q{
                return '/view/genome/druggable-gene/gene-name-report/set/status.html?name=' . ($subject->members)[0]->name;
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
            calculate => q{
                return ($subject->members)[0]->name;
            },
        },
        title => {
            is => 'Text',
            calculate_from => ['subject'],
            calculate => q{
                return ($subject->members)[0]->name . ' druggablegene'
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
                {
                    name => 'name',
                    position => 'content',
                },
                {
                    name => 'alternate_names',
                    position => 'content',
                },
            ],
        },
    ],
};

1;
