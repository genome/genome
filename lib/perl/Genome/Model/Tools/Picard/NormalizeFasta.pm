package Genome::Model::Tools::Picard::NormalizeFasta;

use strict;
use warnings;

use Genome;
use File::Basename;


class Genome::Model::Tools::Picard::NormalizeFasta {
    is  => 'Genome::Model::Tools::Picard',
    has_input => [
        input  => {
            is  => 'String',
            doc => 'The input fasta file to normalize',
        },
        output => {
            is  => 'String',
            doc => 'The output normalized fasta file',
        },
       line_length => {
           is => 'Integer',
           doc => 'The line length to be used for the output fasta file. Default value: 100.',
           is_optional => 1,
       },
       truncate_sequence_names_at_whitespace => {
           is => 'String',
           doc => 'Truncate sequence names at first whitespace. Default value: false.',
           is_optional => 1,
           valid_values => ['true','false'],
       },
    ],
};

sub help_brief {
    'Takes any file that conforms to the fasta format and normalizes it so that all lines of sequence except the last line per named sequence are of the same length.';
}

sub help_detail {
    return <<EOS
    For Picard documentation of this command see:
    http://picard.sourceforge.net/command-line-overview.shtml#NormalizeFasta
EOS
}

sub execute {
    my $self = shift;

    my $cmd = $self->picard_path .'/NormalizeFasta.jar net.sf.picard.reference.NormalizeFasta OUTPUT='. $self->output  .' INPUT='. $self->input;
    if (defined($self->line_length)) {
        $cmd .= ' LINE_LENGTH='. $self->line_length;
    }
    if (defined($self->truncate_sequence_names_at_whitespace)) {
        $cmd .= ' TRUNCATE_SEQUENCE_NAMES_AT_WHITESPACE='. $self->truncate_sequence_names_at_whitespace;
    }
    $self->run_java_vm(
        cmd          => $cmd,
        input_files  => [$self->input],
        output_files => [$self->output],
        skip_if_output_is_present => 0,
    );
    return 1;
}


1;
