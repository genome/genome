package Genome::Model::Tools::Fastq::FilterHomopolymer;

use strict;
use warnings;

use Genome;

my @nucleotides = qw/A T/;
class Genome::Model::Tools::Fastq::FilterHomopolymer {
    is => 'Command',
    has => [
            input => {
                      is => 'Text',
                      doc => 'the input fastq file',
                  },
            output => {
                       is => 'Text',
                       doc => 'the file name to use for the output file',
                   },
            nucleotides => {
                            is => 'List',
                            doc => 'the nucleotide homopolymer runs to remove(default_value='. join(",",@nucleotides) .')',
                            default_value => \@nucleotides,
                        },
            length => {
                       is => 'Number',
                       doc => 'the length of the nucleotide homopolymer run to remove(default_value=10)',
                       default_value => '10',
                  }
    ],
    has_optional => [
                     removed => {
                                 is => 'Text',
                                 doc => 'a fastq file of removed reads',
                             },
    ],
    doc => 'Remove homopolymer reads from a fastq, or pair of paired-end fastqs.',
};

sub help_synopsis {
    return <<EOS
gmt fastq filter-homopolymer --input=lane1.fastq --output=lane1.filtered.fastq

EOS
}

sub help_detail {
    return <<EOS 
Filters a fastq or fastq pair, removing homopolymer reads.
EOS
}

sub execute {
    my $self = shift;

    my $input_fh = IO::File->new($self->input);
    unless ($input_fh) {
        $self->error_message("Failed to open input file " . $self->input . ": $!");
        return;
    }

    my $output_fh = IO::File->new('>'.$self->output);
    unless ($output_fh) {
        $self->error_message("Failed to open output file " . $self->output . ": $!");
        return;
    }
    my $removed_fh;
    if ($self->removed) {
        $removed_fh = IO::File->new('>'.$self->removed);
        unless ($removed_fh) {
            $self->error_message("Failed to open output file " . $self->removed . ": $!");
            return;
        }
    }

    # NOTE: this version does not use ANY of the optional arguments above.

    GETLINE: while (my $header = $input_fh->getline) {
        my $seq = $input_fh->getline;
        my $sep = $input_fh->getline;
        my $qual = $input_fh->getline;
        for my $nucleotide (@{$self->nucleotides}) {
            my $regex = $nucleotide .'{'. $self->length .'}';
            if ($seq =~ /$regex/i) {
                if ($removed_fh) {
                    $removed_fh->print($header,$seq,$sep,$qual);
                }
                next GETLINE;
            }
        }
        $output_fh->print($header,$seq,$sep,$qual);
    }
    $input_fh->close;
    $output_fh->close;
    if ($removed_fh) { $removed_fh->close; }

    return 1;
}

1;

