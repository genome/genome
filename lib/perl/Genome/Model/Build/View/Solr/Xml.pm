package Genome::Model::Build::View::Solr::Xml;

use strict;
use warnings;

use Genome;

class Genome::Model::Build::View::Solr::Xml {
    is => 'Genome::View::Solr::Xml',
    has_field => [
        type => {
            is => 'Text',
            default => 'build'
        },
        display_type => {
            is  => 'Text',
            default => 'Build',
        },
        display_icon_url => {
            is  => 'Text',
            default => 'genome_model_build_32',
        },
        display_url0 => {
            is => 'Text',
            calculate_from => ['subject'],
            calculate => sub { return join ('?id=', '/view/genome/model/build/status.html',$_[0]->id()); },
        },
        display_label1 => {
            is  => 'Text',
            default => 'model',
        },
        display_url1 => {
            is  => 'Text',
            calculate_from => ['subject'],
            calculate => sub { sprintf('/view/genome/model/status.html?id=%s', $_[0]->model->id); },
        },
        display_label2 => {
            is  => 'Text',
            default => 'data directory',
        },
        display_url2 => {
            is  => 'Text',
            calculate_from => ['subject'],
            calculate => sub { sprintf('https://gscweb.gsc.wustl.edu/%s', $_[0]->data_directory); },
        },
        default_aspects => {
            is => 'ARRAY',
            default => [
                {
                    name => 'build_id',
                    position => 'content',
                },
                {
                    name => '__display_name__',
                    position => 'display_title',
                },
                {
                    name => 'date_scheduled',
                    position => 'timestamp',
                },
            ],
        },
    ],
};
