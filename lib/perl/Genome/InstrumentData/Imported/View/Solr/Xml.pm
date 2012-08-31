package Genome::InstrumentData::Imported::View::Solr::Xml;

use strict;
use warnings;

use Genome;

class Genome::InstrumentData::Imported::View::Solr::Xml {
    is => 'Genome::View::Solr::Xml',
    has_field => [
        type => {
            is => 'Text',
            default => 'imported_instrument_data'
        },
        display_type => {
            is  => 'Text',
            default => 'Imported',
        },
        display_icon_url => {
            is  => 'Text',
            default => 'genome_instrumentdata_imported_32',
        },
        display_url0 => {
            is => 'Text',
            calculate_from => ['subject'],
            calculate => sub { return join ('=', '/view/genome/instrument-data/imported/status.html?id',$_[0]->id()); }
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
                    name => 'instrument_data_id',
                    position => 'title',
                },
                {
                    name => 'instrument_data_id',
                    position => 'content',
                },
                {
                    name => 'barcode',
                    position => 'content',
                },
                {
                    name => 'short_run_name',
                    position => 'content',
                },
                {
                    name => 'sra_accession',
                    position => 'content',
                },
                {
                    name => 'sra_sample_id',
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
