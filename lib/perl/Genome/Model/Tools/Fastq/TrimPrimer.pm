package Genome::Model::Tools::Fastq::TrimPrimer;

use strict;
use warnings;

use Genome;
use Bio::Seq;

class Genome::Model::Tools::Fastq::TrimPrimer {
    is => 'Genome::Model::Tools::Fastq',
    has_input => [
            primer_sequence=> {
                       is => 'String',
                       doc => 'the primer sequence to trim',
                       is_optional => 1,
                       default_value => 'GTTTCCCAGTCACGATA',
                   },
            output => {
                       is => 'Text',
                       is_optional => 1,
                       doc => 'the file name to use for the output file',
                   },
            trim_count => {
                       is => 'Number',
                       is_optional => 1,
                       is_output => 1,
                       doc => 'the # of sequences trimmed',
                   },
    ],
    doc => "Remove primer sequence from fastq sequence ends.",
};

sub help_synopsis {
    return <<EOS
gmt fastq trim-primer --fastq_file=lane1.fastq --output=lane1.trimmed.fastq

EOS
}

sub help_detail {
    return <<EOS 
Trims fastq reads from the end.
EOS
}

sub create 
{
    my $class = shift;

    my $self = $class->SUPER::create(@_);
    return $self;
}

sub execute {
    my $self = shift;

    my $input_fh = IO::File->new($self->fastq_file);
    unless ($input_fh) {
        $self->error_message("Failed to open input file " . $self->fastq_file . ": $!");
        return;
    }

    my $output = ($self->output ? $self->output : ($self->fastq_file=~m/(.*.)(\..*)/ ? "$1.TRIMMED$2" : $self->fastq_file . "TRIMMED.txt"));
    my $output_fh = IO::File->new('>'.$output);
    unless ($output_fh) 
    {
        $self->error_message("Failed to open output file " . $output . ": $!");
        return;
    }

    my $screened = "$output.SCREENED";
    my $screened_fh = IO::File->new('>'.$screened);
    unless ($screened_fh) 
    {
        $self->error_message("Failed to open output file " . $screened . ": $!");
        return;
    }

    my $primer_sequence = $self->primer_sequence;
    my $seqobj = Bio::Seq->new(-seq => $primer_sequence);
    my $reverse_sequence = $seqobj->revcom->seq;
    my $trim_count = 0;

    while (my $header = $input_fh->getline) 
    {
        my $seq = $input_fh->getline;
        my $sep = $input_fh->getline;
        my $qual = $input_fh->getline;
        chomp($seq);
        chomp($qual);

       if ($seq=~m/(^$primer_sequence)(.*)/ or $seq=~m/(^$reverse_sequence)(.*)/)
       {
           $screened_fh->print("$seq\n"); #log read before trimming
           $seq = $2 unless (length($2) < 0 and die("The primer sequence is longer than read '$seq'")); 
           $qual = substr($qual, length($1));
           $trim_count++;
       } 

        # TODO: maybe the read names should containt the positions of the bases like _10-50
        $output_fh->print($header,$seq."\n",$sep,$qual."\n");
    }
    $input_fh->close;
    $output_fh->close;
    $screened_fh->close;
    
    $self->trim_count($trim_count);
    return 1;
}

1;

