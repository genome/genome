package Genome::Model::Tools::BedTools::Subtract;

use strict;
use warnings;

use Genome;

class Genome::Model::Tools::BedTools::Subtract {
    is => 'Genome::Model::Tools::BedTools',
    has_input => [
        input_file_a => {
            is => 'Text',
            doc => 'The file of lines with which to subtract b from.',
            shell_args_position => 1,
        },
        input_file_b => {
            is => 'Text',
            doc => 'The file of lines to subtract from a.',
            shell_args_position => 2,
        },
        strandedness => {
            is => 'Text',
            is_optional => 1,
            doc => 'Force strandedness.  That is, only include hits in A that overlap B on the same or opposite strand.  By default, overlaps are subtracted without respect to strand.',
            valid_values => ['same','opposite'],
        },
        minimum_overlap => {
            is => 'Number',
            is_optional => 1,
            doc => 'Minimum overlap required as a fraction of A.'
        },
        output_file => {
            is => 'Text',
            doc => 'The output file to write subtracted output',
        },
    ],
};

sub help_brief {
    "Removes the portion(s) of an interval that is overlapped by another feature(s).",
}

sub help_synopsis {
    my $self = shift;
    return <<"EOS"
  gmt bed-tools subtract --input-file-a genes.bed --input-file-b introns.bed --output-file genes_minus_introns.bed
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
    if ($self->strandedness) {
        if ($self->strandedness eq 'same') {
            $options .= ' -s';
        } elsif ($self->strandedness eq 'opposite') {
            $options .= ' -S';
        }
    }
    if (defined($self->minimum_overlap)) {
        $options .= ' -f '. $self->minimum_overlap;
    }
    my $cmd = $self->bedtools_path .'/bin/subtractBed '. $options .' -a '. $self->input_file_a .' -b '. $self->input_file_b .' > '. $self->output_file;
    Genome::Sys->shellcmd(
        cmd => $cmd,
        input_files => [$self->input_file_a,$self->input_file_b],
        output_files => [$self->output_file],
	skip_if_output_is_present => 0,
    );
    return 1;
}

1;
