package Genome::ModelGroup::View::Solr::Xml;

use strict;
use warnings;

use Genome;

class Genome::ModelGroup::View::Solr::Xml {
    is => 'Genome::View::Solr::Xml',
    has_field => [
        type => {
            is => 'Text',
            default => 'modelgroup'
        },
        display_type => {
            is  => 'Text',
            default => 'Model Group',
        },
        display_icon_url => {
            is  => 'Text',
            default => 'genome_modelgroup_32',
        },
        display_url0 => {
            is => 'Text',
            calculate_from => ['subject'],
            calculate => sub { 
                return join ('?id=', '/view/genome/model-group/status.html',$_[0]->id()); 
            },
        },
        default_aspects => {
            is => 'ARRAY',
            default => [
                {
                    name => '_model_content_for_search',
                    position => 'content',
                },
                {
                    name => '__display_name__',
                    position => 'display_title',
                },
            ]
        }
    ]
};






1;
