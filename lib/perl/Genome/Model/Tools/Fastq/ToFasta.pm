package Genome::Model::Tools::Fastq::ToFasta;

use strict;
use warnings;

use Genome;

class Genome::Model::Tools::Fastq::ToFasta {
    is => 'Genome::Model::Tools::Fastq',
    has_input => [
        fasta_file => {
            is => 'Text',
            doc => 'the output fasta format sequence file',
        },
        quality_file => {
            is => 'Text',
            is_optional => 1,
            is_deprecated => 1,
            doc => 'the output fasta format quality file.  Deprecated since quality conversion not implemented.',
        },
    ],
};

sub execute {
    my $self = shift;
    unless ($self->quality_file) {
        $self->fastQ2A;
    } else {
        $self->bio_convert;
    }
    return 1;
}

sub fastQ2A {
    my $self = shift;

    my $in_fh = Genome::Sys->open_file_for_reading($self->fastq_file);
    my $out_fh = Genome::Sys->open_file_for_writing($self->fasta_file);
    while (my $desc_line = $in_fh->getline) {
        # @HWI-EAS75:1:2:0:345#0/1
        # NCCGCGAGATCGGAAGAGCGGTTCAGCAGGAATGC
        # +HWI-EAS75:1:2:0:345#0/1
        # ENUUUXVUVUUTPTXTSTTQQQVVTTPQVVVVUTQ
        $desc_line =~ s/^\@/\>/;
        my $seq_line = $in_fh->getline;
        my $opt_line = $in_fh->getline;
        my $qual_line = $in_fh->getline;
        print $out_fh $desc_line;
        print $out_fh $seq_line;
    }
    $out_fh->close;
    $in_fh->close;
    return 1;
}

sub bio_convert {
    my $self = shift;
    # TODO: may need a flag to perform quality conversion as well?
    $self->error_message('For quality too, please implement bio_convert method in '. __PACKAGE__);
    die $self->error_message;
}
