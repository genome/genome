package Genome::Model::Tools::CopyNumber::BamWindow;

##############################################################################
#	AUTHOR:		Chris Miller (cmiller@genome.wustl.edu)
#	CREATED:	11/30/2011 by CAM.
#	NOTES:
##############################################################################

use strict;
use Genome;
use IO::File;
use warnings;
use FileHandle;

class Genome::Model::Tools::CopyNumber::BamWindow {
    is => 'Command',
    has => [
    bam_file => {
        is => 'String',
        is_optional => 0,
        doc => 'bam file to count reads from',
    },

    per_library => {
        is => 'Boolean',
        is_optional => 1,
        default => 1,
        doc => 'do counts on a per-library basis',
    },

    output_file => {
        is => 'String',
        is_optional => 0,
        doc => 'output file',
    },

    minimum_mapping_quality => {
        is => 'Integer',
        is_optional => 1,
        doc => 'minimum mapping quality required for a read to be included',
        default => 1,
    },

    extra_params => {
        is => 'String',
        is_optional => 1,
        doc => 'extra parameters to pass to bam-window',
        default => "-s",
    },

    window_size => {
        is => 'Integer',
        is_optional => 1,
        doc => 'size of the bins to count reads in',
        default => 10000,
    },

    per_read_length => {
        is => 'Boolean',
        is_optional => 1,
        doc => 'split different read lengths out into columns',
        default => 1,
    },

    version => {
        is => 'Float',
        is_optional => 1,
        doc => 'Version of bam-window to use',
        default => "0.5",
    },

    ]
};

sub help_brief {
    "Wrapper around bam-window"
}

sub help_detail {
    "Wrapper around bam-window"
}

#########################################################################

sub execute {
    my $self = shift;

    #create the bam-window command
    my $cmd = "/usr/bin/bam-window"; 
    if($self->version eq "0.4"){
        $cmd = $cmd . $self->version;
    }

    $cmd .= " -q " . $self->minimum_mapping_quality;
    $cmd .= " -w " . $self->window_size;

    if($self->per_library){
        $cmd .= " -l";
    }

    if($self->per_read_length){
        $cmd .= " -r"
    }

    $cmd .= " " . $self->extra_params;
    $cmd .= " " . $self->bam_file;
    $cmd .= " >" . $self->output_file;

    my $return = Genome::Sys->shellcmd(
        cmd => "$cmd",
    );
    unless($return) {
        $self->error_message("Failed to execute: Returned $return");
        die $self->error_message;
    }
    return $return;
}
