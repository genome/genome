package Genome::Wiki::Document::View::Solr::Xml;

use strict;
use warnings;

use Genome;

class Genome::Wiki::Document::View::Solr::Xml {
    is => 'Genome::View::Solr::Xml',
    has_field => [
        type => {
            is => 'Text',
            default => 'wiki-page'
        },
        display_type => {
            is => 'Text',
            default => 'Wiki Page', 
        },
        display_icon_url => {
            is  => 'Text',
            default => 'genome_wiki_document_32',
        },
        display_url0 => {
            is => 'Text',
            calculate_from => ['subject'],
            calculate => sub { 

                my $title = $_[0]->title();
                $title =~ s/\%22//g;
                return join ('', $ENV{GENOME_SYS_SERVICES_WIKI_URL}, $title);
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
        default_aspects => {
            is => 'ARRAY',
            default => [
                {
                    name => 'timestamp',
                    position => 'timestamp',
                },
                {
                    name => 'title',
                    position => 'content',
                },
                {
                    name => 'content',
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
