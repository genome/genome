package Genome::Model::Tools::CopyNumber::GraphBamToCnaRegion;

use strict;
use Genome;

class Genome::Model::Tools::CopyNumber::GraphBamToCnaRegion {
    is => 'Command',
    has => [
    cna_file => {
        is => 'String',
        is_optional => 0,
        doc => 'gmt somatic bam-to-cna output file.',
    },
    output_file => {
        is => 'String',
        is_optional => 0,
        doc => 'File to write image to. Should include .png extension.',
    },
    title => {
        is => 'String',
        is_optional => 1,
        doc => 'Title of the graph.',
    },
    color_roi => {
        is => 'Boolean',
        is_optional => 1,
        default => 1,
        doc => 'Whether or not to color ROI points in red',
    },
    chromosome => {
        is => 'String',
        is_optional => 0,
        doc => 'Chromosome of the region to graph',
    },
    roi_start => {
        is => 'Integer',
        is_optional => 0,
        doc => 'The start position of the region of interest on the chromosome.',
    },
    roi_end => {
        is => 'Integer',
        is_optional => 0,
        doc => 'The end position of the region of interest on the chromosome.',
    },
    start => {
        is => 'Integer',
        is_optional => 1,
        doc => 'The start position of the region to be graphed on the chromosome.',
    },	
    end => {
        is => 'Integer',
        is_optional => 1,
        doc => 'The end position of the region to be graphed on the chromosome.',
    },	
    ]
};

sub help_brief {
    "Graphs a region from a BamToCna file"
}

sub help_detail {
    "This script graphs a small region from a BamToCna file. It attempts to output similar to gmt copy-number graph, but without annotation and using the normalization from the BamToCna output"
}

sub execute {
    my $self = shift;
    my $rlibrary = "BamToCnaGraph.R";

    #TODO some sanity checks on the input
    
    my $color = $self->color_roi ? "T" : "F";

    #Call R for graphing 
    my $graph_cmd = sprintf("graph_cna_region(cna_file='%s',output_file='%s',chromosome='%s',region_start=%d,region_end=%d,roi_start=%d,roi_end=%d,color_roi=%s);",$self->cna_file,$self->output_file,$self->chromosome,$self->start,$self->end,$self->roi_start,$self->roi_end,$color);
    my $graph_rcall = Genome::Model::Tools::R::CallR->create(command=>$graph_cmd,library=>$rlibrary);
    $graph_rcall->execute;

}

1;
