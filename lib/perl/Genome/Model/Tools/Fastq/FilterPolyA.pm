package Genome::Model::Tools::Fastq::FilterPolyA;

use strict;
use warnings;

use Genome;
use Bio::SeqIO;
use Bio::Seq::Quality;

class Genome::Model::Tools::Fastq::FilterPolyA {
    is => 'Command',
    has => [    
            input1 => {
                is => 'Text',
                doc => 'the input fastq file, or the first of a pair', 
            },
            output1 => {
                is => 'Text',
                doc => 'the file name to use for the output file, or the first of a pair',
            },
    ],
    has_optional => [
            input2 => {
                is => 'Text',
                doc => 'the second input fastq file when filtering paired-end data',
            },
            output2 => {
                is => 'Text',
                doc => 'the second output fastq file when filtering paired-end data',
            },
            removed1 => {
                is => 'Text',
                is_optional => 1,
                doc => 'write all removed reads from input1 into anther fastq',
            },            
            removed2 => {
                is => 'Text',
                is_optional => 1,
                doc => 'write all removed reads from input2 into anther fastq',
            },            
    ],
    doc => 'Remove poly-a reads from a fastq, or pair of paired-end fastqs.',
};

sub help_synopsis {
    return <<EOS
gmt fastq filter-poly-a --input1=lane1.fastq --output1=lane1.filtered.fastq

gmt fastq filter-poly-a --input1=lane1-fwd.fastq --input2=lane2-rev.fastq --output1=lane1-fwd.filtered.fastq --output2=lane1-rev.filtered.fastq
(not implemented)
EOS
}

sub help_detail {
    return <<EOS 
Filters a fastq or fastq pair, removing poly-A reads.
EOS
}

sub execute {
    my $self = shift;

    my $input1_fh = IO::File->new($self->input1);
    unless ($input1_fh) {
        $self->error_message("Failed to open input file " . $self->input1 . ": $!");
        return;
    }

    my $output1_fh = IO::File->new('>'.$self->output1);
    unless ($output1_fh) {
        $self->error_message("Failed to open output file " . $self->output1 . ": $!");
        return;
    }

    # NOTE: this version does not use ANY of the optional arguments above.

    while (my $header = $input1_fh->getline) {
        my $seq = $input1_fh->getline;
        my $sep = $input1_fh->getline;
        my $qual = $input1_fh->getline;
        if($seq =~ /A{15}/i) { 
            # 15 As in a row somewhere in the read: skip
            next;
        }
        $output1_fh->print($header,$seq,$sep,$qual);
    }

    $input1_fh->close;
    $output1_fh->close;

    return 1;
}

1;

