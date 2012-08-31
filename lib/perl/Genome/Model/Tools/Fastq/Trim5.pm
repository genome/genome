package Genome::Model::Tools::Fastq::Trim5;

use strict;
use warnings;

use Genome;

class Genome::Model::Tools::Fastq::Trim5 {
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
            length => {
                       is => 'Number',
                       doc => 'the length of the nucleotide homopolymer run to remove(default_value=9)',
                       default_value => '9',
                  }
    ],
    has_optional => [
                     removed => {
                                 is => 'Text',
                                 doc => 'a fastq file with the trimmed or excised sequence',
                             },
                 ],
    doc => "Remove 5' base pairs from fastq sequence.",
};

sub help_synopsis {
    return <<EOS
gmt fastq trim5 --input=lane1.fastq --output=lane1.trimmed.fastq

EOS
}

sub help_detail {
    return <<EOS 
Trims fastq reads from the 5' end.
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

    # NOTE: this version does not use ANY of the optional arguments above.
    my $removed_fh;
    if ($self->removed) {
        $removed_fh = IO::File->new('>'.$self->removed);
        unless ($removed_fh) {
            $self->error_message("Failed to open removed file " . $self->removed . ": $!");
            return;
        }
    }

    while (my $header = $input_fh->getline) {
        my $seq = $input_fh->getline;
        my $sep = $input_fh->getline;
        my $qual = $input_fh->getline;
        chomp($seq);
        chomp($qual);
        my @nucleotides = split("",$seq);
        my @quality_values = split("",$qual);

        #This is for 5' only
        # TODO: maybe the read names should containt the positions of the bases like _10-50
        my @spliced_nucleotides = splice(@nucleotides,0,$self->length);
        my @spliced_quality_values = splice(@quality_values,0,$self->length);
        $output_fh->print($header,join("",@nucleotides)."\n",$sep,join("",@quality_values)."\n");
        if ($removed_fh) {
            # TODO: maybe the read names should containt the positions of the bases like _1-9
            chomp($header);
            chomp($seq);
            $removed_fh->print($header ."-5'-". $self->length ."\n",
                               join("",@spliced_nucleotides)."\n",
                               $sep ."-5'-". $self->length ."\n",
                               join("",@spliced_quality_values)."\n");
        }
    }
    $input_fh->close;
    $output_fh->close;
    if ($removed_fh) { $removed_fh->close; }
    return 1;
}

1;

