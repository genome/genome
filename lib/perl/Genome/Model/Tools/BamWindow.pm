package Genome::Model::Tools::BamWindow;

use strict;
use warnings;

use Genome;

class Genome::Model::Tools::BamWindow {
    is => 'Command::V2',
    has_input => {
        version => {
            is => 'Text',
            default => '0.4',
            doc => 'Version of bam-window to run',
        },
        options => {
            is => 'Text',
            doc => 'String of command line options to pass to bam-window.  Superseeds all other options',
        },
        bam_file => {
            is => 'Path',
            doc => 'Bam file to window',
        },
        output_file => {
            is => 'Path',
            is_output => 1,
            doc => 'Window output file',
        },
        quality => {
            is => 'Number',
            doc => 'filtering reads with mapping quality less than INT [0]',
            is_optional => 1,
        },
        window_size => {
            is => 'Number',
            doc => 'window size to count reads within [1000]',
            is_optional => 1,
        },
        paired_reads_only => {
            is => 'Boolean',
            doc => 'only include paired reads',
            is_optional => 1,
        },
        properly_paired_reads_only => {
            is => 'Boolean',
            doc => 'only include properly paired reads',
            is_optional => 1,
        },
        leftmost_only => {
            is => 'Boolean',
            doc => 'only count a read as in the window if its leftmost mapping position is within the window',
            is_optional => 1,
        },
        per_library => {
            is => 'Boolean',
            doc => 'output a column for each library in each window',
            is_optional => 1,
        },
        per_read_length => {
            is => 'Boolean',
            doc => 'output a column for each read length in each window',
            is_optional => 1,
        },
        probability => {
            is => 'Number',
            doc => 'probability of reporting a read [1.000000]',
            is_optional => 1,
        },
    },
};

sub execute {
    my $self = shift;
    
    my $base_cmd = "bam-window"; #TODO: choose this based on $self->version
    my $bam_file = $self->bam_file;

    my $options_string = $self->options;
    unless($options_string){
        $self->_get_options_string;
    }
    my $output_file = $self->output_file;

    my $cmd = join(" ", $base_cmd, $bam_file, $options_string, " > $output_file");
    system($cmd);
    if($@){
        $self->error_message("Error running command $cmd: $@");
        return 0;
    }

    return 1;
}

sub _get_options_string {
    my $self = shift;
    my $options_string = "";

    if ($self->quality){
        $options_string = join(" ", $options_string, "-q", $self->quality);
    }

    if($self->window_size){
        $options_string = join(" ", $options_string, "-w", $self->window_size);
    }

    if($self->paired_reads_only){
        $options_string = join(" ", $options_string, "-p");
    }

    if($self->properly_paired_reads_only){
        $options_string = join(" ", $options_string, "-P");
    }

    if($self->leftmost_only){
        $options_string = join(" ", $options_string, "-s");
    }

    if($self->per_library){
        $options_string = join(" ", $options_string, "-l");
    }

    if($self->per_read_length){
        $options_string = join(" ", $options_string, "-r");
    }

    if($self->probability){
        $options_string = join(" ", $options_string, "-d", $self->probability);
    }

    return $options_string;
}

1;
