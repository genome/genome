package Genome::Model::Tools::Fastq::Trimmer;

use strict;
use warnings;

use Genome;

class Genome::Model::Tools::Fastq::Trimmer {
    is => 'Command',
    has => [
            input_fastq => {
                      is => 'Text',
                      doc => 'the input fastq file',
                  },
            output_fastq => {
                       is => 'Text',
                       doc => 'the file name to use for the output file',
                   },
            five_prime => {
                is => 'Number',
                doc => 'the number of base pair to trim from the 5\'-prime end of the read',
                default_value => 0,
                is_optional => 1,
            },
            three_prime => {
                is => 'Number',
                doc => 'the number of base pair to trim from the 3\'-prime end of the read',
                default_value => 0,
                is_optional => 1,
            },
    ],
    doc => "Remove base pairs from fastq sequence ends.",
};

sub help_synopsis {
    return <<EOS
gmt fastq trimmer --input=lane1.fastq --output=lane1.trimmed.fastq --five=10 --three=25

EOS
}

sub help_detail {
    return <<EOS 
Trims fastq reads from the end.
EOS
}

sub execute {
    my $self = shift;

    my $input_fh = IO::File->new($self->input_fastq);
    unless ($input_fh) {
        $self->error_message("Failed to open input fastq file " . $self->input_fastq . ": $!");
        return;
    }

    my $output_fh = IO::File->new('>'.$self->output_fastq);
    unless ($output_fh) {
        $self->error_message("Failed to open output file " . $self->output_fastq . ": $!");
        return;
    }

    my $min_seq_length = $self->five_prime + $self->three_prime;
    while (my $header = $input_fh->getline) {
        my $seq = $input_fh->getline;
        my $sep = $input_fh->getline;
        my $qual = $input_fh->getline;
        chomp($seq);
        chomp($qual);

        my $seq_length = length($seq);
        if ($seq_length <= $min_seq_length) {
            die('fastq sequence for '. $header .' is only '. $seq_length .' base pair and min read length calculated as '. $min_seq_length);
        }
        my $offset = 0;
        if ($self->five_prime) {
            $offset = $self->five_prime;
        }
        my $length = $seq_length;
        if ($self->three_prime) {
            $length = $seq_length - $self->three_prime - $self->five_prime;
        }
        $seq = substr($seq,$offset,$length);
        $qual = substr($qual,$offset,$length);
        $output_fh->print($header . $seq ."\n". $sep . $qual ."\n");
    }
    $input_fh->close;
    $output_fh->close;
    return 1;
}

1;

