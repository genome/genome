package Genome::Taxon::View::Solr::Xml;

use strict;
use warnings;

use Genome;

class Genome::Taxon::View::Solr::Xml {
    is => 'Genome::View::Solr::Xml',
    has_field => [
        type => {
            is => 'Text',
            default => 'taxon'
        },
        display_type => {
            is  => 'Text',
            default => 'Taxon',
        },
        display_icon_url => {
            is  => 'Text',
            default => 'genome_taxon_32',
        },
        display_url0 => {
            is => 'Text',
            calculate_from => ['subject'],
            calculate => sub { return join ('?id=', '/view/genome/taxon/status.html',$_[0]->id()); },
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
        default_aspects => {
            is => 'ARRAY',
            default => [
                {
                    name => 'domain',
                    position => 'content',
                },
                {
                    name => 'species_name',
                    position => 'content',
                },
                {
                    name => 'strain_name',
                    position => 'content',
                },
                {
                    name => 'species_latin_name',
                    position => 'content',
                },
                {
                    name => 'ncbi_taxon_id',
                    position => 'content',
                },
                {
                    name => '__display_name__',
                    position => 'display_title',
                },
            ],
        }
    ]
};

1;
