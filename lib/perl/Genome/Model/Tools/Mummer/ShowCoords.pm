package Genome::Model::Tools::Mummer::ShowCoords;

use strict;
use warnings;

use Genome;

class Genome::Model::Tools::Mummer::ShowCoords {
    is => 'Genome::Model::Tools::Mummer',
    has => [
        input_delta_file => {
            is => 'Text',
            doc => 'Alignment file in delta format',
        },
        output_file => {
            is => 'Text',
            doc => 'File to direct output to',
        },
        b => {
            is => 'Boolean',
            doc => 'Merges overlapping alignments regardless of match dir or frame and does not display any idenitity information',
            is_optional => 1,
        },
        c => {
            is => 'Boolean',
            doc => 'Include percent coverage information in the output',
            is_optional => 1,
        },
        l => {
            is => 'Boolean',
            doc => 'Include the sequence length information in the output',
            is_optional => 1,
        },
        o => {
            is => 'Boolean',
            doc => 'Annotate maximal alignments between two sequences',
            is_optional => 1,
        },
       'q' => {
            is => 'Boolean',
            doc => 'Sort output lines by query IDs and coordinates',
            is_optional => 1,
        },
        r => {
            is => 'Boolean',
            doc => 'Sort output lines by reference IDs and coordinates',
            is_optional => 1,
        },
        B => {
            is => 'Boolean',
            doc => 'Switch output to btab format',
            is_optional => 1,
        },
        H => {
            is => 'Boolean',
            doc => 'Do not print the output header',
            is_optional => 1,
        },
        I => {
            is => 'float',
            doc => 'Set minimum percent identity to display',
            is_optional => 1,
        },
        M => {
            is => 'Integer',
            doc => 'Set minimum alignment length to display',
            is_optional => 1,
        },
        T => {
            is => 'Boolean',
            doc => 'Switch output to tab-delimited format',
            is_optional => 1,
        },
    ],
};

sub help_brief {
    'Tool to run show-coords on nucmer generated alignments in delta file format';
}

sub help_detail {
    return <<EOS
gmt mummer nucer show-coords --input-delta-file OUT.delta --output-file alignments.txt -r -c -l -T
EOS
}

sub execute {
    my $self = shift;

    if ( not -s $self->input_delta_file ) {
        $self->error_message("Failed to find input delta file: ".$self->input_delta_file);
        return;
    }

    if ( -e $self->output_file ) {
        $self->debug_message('Removing existing output file: '.$self->output_file);
        unlink $self->output_file;
    }

    my $output_dir = File::Basename::dirname($self->output_file);
    if ( $output_dir and not -d $output_dir ) {
        $self->error_message("Failed to find output dir: $output_dir");
        return;
    }

    my $versioned_show_coords = $self->path_for_version.'/show-coords';
    if ( not -s $versioned_show_coords ) {
        $self->error_message('Failed to find show-coords in version: '.$self->path_for_version);
        return;
    }

    my $params;
    $params .= $self->stringify_boolean_params;
    $params .= ' -I '.$self->I if $self->I;
    $params .= ' -L '.$self->M if $self->M;

    my $cmd = $versioned_show_coords.$params.' '.$self->input_delta_file.' > '.$self->output_file;

    $self->debug_message("Running show-coords with command: $cmd");
    
    my %run_params = (
        cmd          => $cmd,
        #input_files  => [ $self->input_delta_file ],
        output_files => [ $self->output_file ], 
    );
    my $rv = eval{ Genome::Sys->shellcmd( %run_params ); };
    if ( not $rv ) {
        $self->error_message("Failed to run show-coords with command: $cmd");
        return;
    }

    return 1;
}

sub stringify_boolean_params {
    my $self = shift;

    my $string;
    for my $param ( qw/ b c l o q r B H T / ) {
        if ( defined $self->$param ) {
            $string .= $param;
        }
    }
    return ' -'.$string if $string;
    return '';
}

1;
