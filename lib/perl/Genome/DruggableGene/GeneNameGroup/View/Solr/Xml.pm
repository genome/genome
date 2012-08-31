package Genome::DruggableGene::GeneNameGroup::View::Solr::Xml;

use strict;
use warnings;

use Genome;

class Genome::DruggableGene::GeneNameGroup::View::Solr::Xml {
    is => 'Genome::View::Solr::Xml',
    has_field => [
        type => {
            is => 'Text',
            default => 'gene-name-group'
        },
        display_type => {
            is  => 'Text',
            default => 'GeneNameGroup',
        },
        display_icon_url => {
            is  => 'Text',
            default => 'genome_druggable-gene_gene-name-group_32',
        },
        display_url0 => {
            is => 'Text',
            calculate_from => ['subject'],
            calculate => q{
                return '/view/genome/druggable-gene/gene-name-group/status.html?name=' . $subject->name;
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
                return $subject->name
            },
        },
        title => {
            is => 'Text',
            calculate_from => ['subject'],
            calculate => q{
                return join(' ',map{$_->name}$subject->genes) .
                join(' ', map{$_->alternate_names}$subject->genes) .
                ' druggablegene'
            },
        },
        default_aspects => {
            is => 'ARRAY',
            default => [
                {
                    name => 'name',
                    position => 'content',
                },
            ],
        },
    ],
};

1;
