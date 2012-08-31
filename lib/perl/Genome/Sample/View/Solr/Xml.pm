package Genome::Sample::View::Solr::Xml;

use strict;
use warnings;

use Genome;

class Genome::Sample::View::Solr::Xml {
    is => 'Genome::View::Solr::Xml',
    has_field => [
        type => {
            is => 'Text',
            default => 'sample'
        },
        display_type => {
            is  => 'Text',
            default => 'Sample',
        },
        display_icon_url => {
            is  => 'Text',
            default => 'genome_sample_32',
        },
        display_url0 => {
            is => 'Text',
            calculate_from => ['subject'],
            calculate => sub { return join ('?id=', '/view/genome/sample/status.html',$_[0]->id()); },
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
                    name => 'extraction_label',
                    position => 'content',
                },
                {
                    name => 'extraction_type',
                    position => 'content',
                },
                {
                    name => 'extraction_desc',
                    position => 'content',
                },
                {
                    name => 'individual_common_name',
                    position => 'content',
                },
                {
                    name => 'cell_type',
                    position => 'content',
                },
                {
                    name => 'tissue_label',
                    position => 'content',
                },
                {
                    name => 'tissue_desc',
                    position => 'content',
                },
                {
                    name => 'organ_name',
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
