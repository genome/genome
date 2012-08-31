package Genome::InstrumentData::FlowCell::View::Solr::Xml;

use strict;
use warnings;

use Genome;

class Genome::InstrumentData::FlowCell::View::Solr::Xml {
    is => 'Genome::View::Solr::Xml',
    has_field => [
        type => {
            is => 'Text',
            default => 'illumina_run'
        },
        display_type => {
            is  => 'Text',
            default => 'Flow Cell',
        },
        display_icon_url => {
            is  => 'Text',
            default => 'genome_instrumentdata_flowcell_32',
        },
        display_url0 => {
            is => 'Text',
            calculate_from => ['subject'],
            calculate => sub { return join ('=', '/view/genome/instrument-data/flow-cell/status.html?id',$_[0]->id()); }
        },
        display_label1 => {
            is  => 'Text',
            default => 'production run',
        },
        display_url1 => {
            is  => 'Text',
            calculate_from => ['subject'],
            calculate => sub {
                my $flow_cell_id = $_[0]->flow_cell_id();
                return 'none' if !$flow_cell_id;
                return '/solexa/equipment/flowcell/' . $flow_cell_id;
            },
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
                    name => 'flow_cell_id',
                    position => 'title',
                },
                {
                    name => 'instrument_data_ids',
                    position => 'content',
                },
                {
                    name => 'flow_cell_id',
                    position => 'content',
                },
                {
                    name => 'machine_name',
                    position => 'content',
                },
                {
                    name => 'run_name',
                    position => 'content',
                },
                {
                    name => 'run_type',
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
