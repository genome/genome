package Genome::Model::Tools::Fasta::SlidingWindows;

use strict;
use warnings;

use Genome;

class Genome::Model::Tools::Fasta::SlidingWindows (
    is => ['Genome::Model::Tools::Fasta'],
    has => [
        output_file => {
            is => 'Text',
            doc => 'The output FASTA file',
        },
    ],
    has_optional => [
        window_size => {
            is => 'Integer',
            default_value => 150,
            doc => 'The size of the window or sub-sequence to create.',
        },
        window_offset => {
            is => 'Integer',
            default_value => 5,
            doc => 'The amount of offset between each sliding window.',
        },
        direction => {
            is => 'Text',
            valid_values => ['fwd','rev','both'],
            default_value => 'both',
            doc => 'The direction to slide window.  fwd = left to right, rev = right to left',
        },
        #TODO Implement reverse complement when using the rev strand.
        #For primer design however this is bad,  the primers will anneale to each other
    ],
);

sub execute {
    my $self = shift;
    my $reader = $self->get_fasta_reader($self->fasta_file);
    unless ($reader) {
        $self->error_message('Failed to open FASTA reader for file '. $self->fasta_file);
        die($self->error_message);
    }
    my $writer = $self->get_fasta_writer($self->output_file);
    unless ($writer) {
        $self->error_message('Failed to open FASTA writer for file '. $self->output_file);
        die($self->error_message);
    }
    while (my $seq = $reader->next_seq) {
        my $length = $seq->length;
        if ($self->direction eq 'fwd' || $self->direction eq 'both') {
            for (my $i = 1; $i <= ((($length - $self->window_size) + 1)); $i += $self->window_offset) {
                my $j = $i + $self->window_size - 1;
                my $substring = $seq->subseq($i,$j);
                my $substring_id = $seq->id .':'. $i .'-'. $j .'_fwd';
                my $substring_seq = Bio::PrimarySeq->new (
                    -seq => $substring,
                    -id  => $substring_id,
                    -alphabet => 'dna',
                );
                $writer->write_seq($substring_seq);
            }
        }
        if ($self->direction eq 'rev' || $self->direction eq 'both') {
            for (my $j = $length; $j >= ($self->window_size); $j -= $self->window_offset) {
                my $i = (($j - $self->window_size) + 1);
                my $substring = $seq->subseq($i,$j);
                my $substring_id = $seq->id .':'. $i .'-'. $j .'_rev';
                my $substring_seq = Bio::PrimarySeq->new (
                    -seq => $substring,
                    -id  => $substring_id,
                    -alphabet => 'dna',
                );
                $writer->write_seq($substring_seq);
            }
        }
    }
    return 1;
}


1;
