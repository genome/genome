package Genome::Project::View::Solr::Xml;

use strict;
use warnings;

use Genome;

class Genome::Project::View::Solr::Xml {
    is => 'Genome::View::Solr::Xml',
    has_field => [
        type => {
            is => 'Text',
            default => 'project'
        },
        display_type => {
            is  => 'Text',
            default => 'Project',
        },
        display_icon_url => {
            is  => 'Text',
            default => 'genome_project_32',
        },
        display_url0 => {
            is => 'Text',
            calculate_from => ['subject'],
            calculate => sub { 
                return join ('?id=', '/view/genome/project/status.html', $_[0]->id()); 
            },
        },
        display_url1 => {},
        display_label1 => {},
        display_url2 => {},
        display_label2 => {},
        display_url3 => {},
        display_label3 => {},
        default_aspects => {
            is => 'ARRAY',
            default => [
                {
                    name => 'creator',
                    position => 'content',
                },
                {
                    name => '__display_name__',
                    position => 'display_title',
                },
                {   name => 'id',
                    position => 'content',
                },
                {
                    name => 'name',
                    position => 'content',
                },
            ],
        }
    ]
};

