package Genome::Model::Tools::BedTools::Map;

use strict;
use warnings;

use Genome;

my $DEFAULT_OPERATION = 'sum';
my $DEFAULT_FORCE_STRANDEDNESS = 0;
my $DEFAULT_COLUMN = 5;

class Genome::Model::Tools::BedTools::Map {
    is => 'Genome::Model::Tools::BedTools',
    has_input => [
        input_file_a => {
            is => 'Text',
            doc => 'The file of lines with which to find overlaps.',
            shell_args_position => 1,
        },
        input_file_b => {
            is => 'Text',
            doc => 'The file in which to find overlaps for each line of input-file-a.',
            shell_args_position => 2,
        },
        column => {
            is => 'Number',
            doc => 'The column from the B file to map onto intervals in A',
            default_value => $DEFAULT_COLUMN,
            is_optional => 1,
        },
        force_strandedness => {
            is => 'Boolean',
            is_optional => 1,
            default_value => $DEFAULT_FORCE_STRANDEDNESS,
            doc => 'Force strandedness.  That is, only include hits in A that overlap B on the same strand.'
        },
        minimum_overlap => {
            is => 'Number',
            is_optional => 1,
            doc => 'Minimum overlap required as a fraction of A.'
        },
        output_file => {
            is => 'Text',
            doc => 'The output file to write intersection output',
        },
        operation => {
            is => 'Text',
            doc => 'Operation that should be applied to -c',
            valid_values => ['sum','min','max','mean','median','collapse','distinct','count','count_distinct'],
            default => $DEFAULT_OPERATION,
        },
        reciprocal => {
            is => 'Boolean',
            doc => 'Require that the fraction overlap be reciprocal for A and B',
            default => 0,
        },
        header => {
            is => 'Boolean',
            doc => 'Print the header from the A file prior to results',
            default => 0,
        },
    ],
};

sub help_brief {
    "Apply a function to a column from B intervals that overlap A.",
}

sub help_synopsis {
    my $self = shift;
    return <<"EOS"
  gmt bed-tools map --input-file-a example.bam --input-file-b other-list.bed --output-file intersection.bam --operation distinct
EOS
}

sub help_detail {                           
    return <<EOS
More information about the BedTools suite of tools can be found at http://code.google.com/p/bedtools/. 
EOS
}

sub execute {
    my $self = shift;
    
    my $options = '';
    if ($self->force_strandedness) {
        $options .= ' -s';
    }
    
    if (defined($self->minimum_overlap)) {
        $options .= ' -f '. $self->minimum_overlap;
        if ($self->reciprocal) {
            $options .= ' -r ';
        }
    }
    if ($self->operation) {
        $options .= ' -o '.$self->operation;
    }
    if ($self->header) {
        $options .= " -header ";
    }
    if ($self->column) {
        $options .= " -c ".$self->column;
    }
    my $cmd = $self->bedtools_path .'/bin/mapBed '. $options .' '. $self->input_file_a .' -b '. $self->input_file_b .' > '. $self->output_file;
    Genome::Sys->shellcmd(
        cmd => $cmd,
        input_files => [$self->input_file_a,$self->input_file_b],
        output_files => [$self->output_file],
        skip_if_output_is_present => 0,
        allow_zero_size_output_files => 1,
    );
    return 1;
}

1;
