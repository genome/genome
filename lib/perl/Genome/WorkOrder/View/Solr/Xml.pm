package Genome::WorkOrder::View::Solr::Xml;

use strict;
use warnings;

use Genome;

class Genome::WorkOrder::View::Solr::Xml {
    is => 'Genome::View::Solr::Xml',
    has_field => [
        type => {
            is => 'Text',
            default => 'work-order'
        },
        default_aspects => {
            is => 'ARRAY',
            default => [
                {
                    name => 'barcode',
                    position => 'content',
                },
                {
                    name => 'name',
                    position => 'content',
                },
                {
                    name => 'pipeline',
                    position => 'content',
                },
                {
                    name => 'sample_description',
                    position => 'content',
                },
                {
                    name => 'project',
                    position => 'content',
                    perspective => 'default',
                    toolkit => 'xml',
                    aspects => [
                        'id',
                        'name',
                    ]
                },
                {
                    name => '__display_name__',
                    position => 'display_title',
                }
            ]
        },
        display_type => {
            is  => 'Text',
            default => 'Work Order',
        },
        display_icon_url => {
            is  => 'Text',
            default => 'genome_workorder_32',
        },
        display_url0 => {
            is => 'Text',
            calculate_from => ['subject'],
            calculate => sub {
                return join ('?id=', '/view/genome/work-order/status.html', $_[0]->id());
            },
        },
        display_label1 => {
            is  => 'Text',
            default => '',
        },
        display_url1 => {
            is  => 'Text',
            calculate_from => ['subject'],
            calculate => sub {
            },
        },
        display_label2 => {
            is  => 'Text',
            default => '',
        },
        display_url2 => {
            is  => 'Text',
            calculate_from => ['subject'],
            calculate => sub {
            },
        },
        display_label3 => {
            is  => 'Text',
            default => '',
        },
        display_url3 => {
            is  => 'Text',
            calculate_from => ['subject'],
            calculate => sub {
            },
        }
]};



1;

