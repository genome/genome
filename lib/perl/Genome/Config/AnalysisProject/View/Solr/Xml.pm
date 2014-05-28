package Genome::Config::AnalysisProject::View::Solr::Xml;

use strict;
use warnings;

use Genome;

class Genome::Config::AnalysisProject::View::Solr::Xml {
    is => 'Genome::View::Solr::Xml',
    has_field => [
        type => {
            is => 'Text',
            value => 'analysis_project',
        },
        display_type => {
            is => 'Text',
            value => 'Analysis Project',
        },
        display_icon_url => {
            is => 'Text',
            value => 'genome_project_32',
        },
        display_url0 => {
            is => 'Text',
            calculate_from => ['subject'],
            calculate => sub {
                join('?id=', '/view/genome/config/analysis-project/status.html',$_[0]->id())
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
            value => [
                {
                    name => 'created_by',
                    position => 'content',
                },
                {
                    name => '__display_name__',
                    position => 'display_title',
                },
                {
                    name => 'id',
                    position => 'content',
                },
                {
                    name => 'name',
                    position => 'content',
                },
                {
                    name => 'created_at',
                    position => 'timestamp',
                },
            ],
        }
    ],
};

1;
