package Genome::Individual::View::Solr::Xml;

use strict;
use warnings;

use Genome;

class Genome::Individual::View::Solr::Xml {
    is => 'Genome::View::Solr::Xml',
    has_field => [
        type => {
            is => 'Text',
            default => 'individual'
        },
        display_type => {
            is  => 'Text',
            default => 'Individual',
        },
        display_icon_url => {
            is  => 'Text',
            default => 'genome_individual_32',
        },
        display_url0 => {
            is => 'Text',
            calculate => q { 
                    my $subject = $self->subject;
                    return join ('?id=', '/view/genome/individual/status.html',$subject->individual_id()); 
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
                    name => 'common_name',
                    position => 'title',
                },
                {
                    name => 'name',
                    position => 'content',
                },
                {
                    name => 'gender',
                    position => 'content',
                },
                {
                    name => 'upn',
                    position => 'content',
                },
                {
                    name => '__display_name__',
                    position => 'display_title',
                },
            ]
        },
    ]
};

1;
