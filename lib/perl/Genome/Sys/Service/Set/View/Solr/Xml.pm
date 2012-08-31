package Genome::Sys::Service::Set::View::Solr::Xml;

use strict;
use warnings;

use Genome;

class Genome::Sys::Service::Set::View::Solr::Xml {
    is => 'Genome::View::Solr::Xml',
    has_field => [
        type => {
            is => 'Text',
            default => 'service'
        },
        display_type => {
            is  => 'Text',
            default => 'Service',
        },
        display_icon_url => {
            is  => 'Text',
            default => '',
        },
        display_url0 => {
            is => 'Text',
            calculate_from => ['subject'],
            calculate => q{
                return '/view/genome/druggable-gene/drug-name-report/set/status.html?name=' . $subject->name();
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
                    name => 'name',
                    position => 'title',
                },
                {
                    name => 'host',
                    position => 'content',
                },
                {
                    name => 'restart_command',
                    position => 'content',
                },
                {
                    name => 'stop_command',
                    position => 'content',
                },
                {
                    name => 'log_path',
                    position => 'content',
                },
                {
                    name => 'pid_name',
                    position => 'content',
                },
                {
                    name => '__display_name__',
                    position => 'display_title',#text for display_url0
                },
            ],
        },
    ],
};

1;
