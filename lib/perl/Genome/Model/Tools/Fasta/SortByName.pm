package Genome::Model::Tools::Fasta::SortByName;

use strict;
use warnings;

use Genome;

class Genome::Model::Tools::Fasta::SortByName {
    is => 'Command::V2',
    has => [
        input_fasta => {
            is => 'FilePath',
            is_input => 1,
            doc => 'Input fasta file',
        },
    ],
    has_optional => [
        sorted_fasta => {
            is => 'FilePath',
            is_input => 1,
            is_output => 1,
            doc => 'Sorted output fasta',
        },
    ],
    doc => 'Sorts a fasta by sequence name',
};

sub help_detail {
    return 'Sorts the sequences in the provided fasta by name and writes result to another file';
}

sub execute {
    my $self = shift;

    unless ($self->sorted_fasta) {
        $self->sorted_fasta($self->input_fasta . '.sorted');
        $self->debug_message("No sorted fasta file provided, defaulting to " . $self->sorted_fasta . "!");
    }

    $self->debug_message("Sorting fasta " . $self->input_fasta . " by sequence name and putting results in " . $self->sorted_fasta);

    my $input_fh = Genome::Sys->open_file_for_reading($self->input_fasta);
    my $output_fh = Genome::Sys->open_file_for_writing($self->sorted_fasta);

    my %seqs;
    my $current_seq;
    while (my $line = $input_fh->getline) {
        chomp $line;
        if ($line =~ /^>/) {
            $current_seq = $line;
        }
        else {
            $seqs{$current_seq} .= $line;
        }
    }

    for my $seq_name (sort keys %seqs) {
        my $sequence = $seqs{$seq_name};
        $output_fh->print($seq_name . "\n");
        for (my $i = 0; $i < (length $sequence); $i += 80) {
            my $sub_sequence = substr($sequence, $i, 80);
            $output_fh->print($sub_sequence . "\n");
        }
    }

    $self->debug_message("Successfully sorted fasta file " . $self->input_fasta . " by sequence name!");
    return 1;
}

1;

