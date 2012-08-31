package Genome::InstrumentData::Solexa::View::Solr::Xml;

use strict;
use warnings;

use Genome;

class Genome::InstrumentData::Solexa::View::Solr::Xml {
    is => 'Genome::View::Solr::Xml',
    has_field => [
        type => {
            is => 'Text',
            default => 'solexa_instrument_data'
        },
        display_type => {
            is  => 'Text',
            default => 'Solexa',
        },
        display_icon_url => {
            is  => 'Text',
            default => 'genome_instrumentdata_32',
        },
        display_url0 => {
            is => 'Text',
            calculate_from => ['subject'],
            calculate => sub { return join ('=', '/view/genome/instrument-data/solexa/status.html?id',$_[0]->id()); }
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
                    name => 'id',
                    position => 'title',
                },
                {
                    name => 'id',
                    position => 'content',
                },
                {
                    name => '__display_name__',
                    position => 'display_title',
                },
            ],
        },
    ],
};

1;
