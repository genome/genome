package Genome::Model::PhenotypeCorrelation::Command::SortVepOutput;

use Genome;
use Sort::Naturally qw/nsort/;

use strict;
use warnings;

class Genome::Model::PhenotypeCorrelation::Command::SortVepOutput {
    is => "Command::V2",
    doc => "Sort Vep output files by chrom/pos",
    has => [
        input_file => {
            is => "File",
            doc => "Path to Vep output file",
            is_input => 1,
        },
        output_file => {
            is => "File",
            doc => "Path to the output file to write.",
            is_input => 1,
            is_output => 1,
        },
    ],
};

sub help_synopsis {
    return <<"EOS"
EOS
}

sub help_detail {
    return <<"EOS"
EOS
}

sub execute {
    my $self = shift;

    my $input = $self->input_file;
    my $output = $self->output_file;

    `(awk '{if (/^#/) {print > "/dev/stderr"} else {print}}' $input | sort -k2,2V -s) > $output 2>&1`;

    return 1;
}

1;
