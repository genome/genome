package Genome::Model::Tools::RefCov::BamHeaderToBed;

use strict;
use warnings;

use Genome;

class Genome::Model::Tools::RefCov::BamHeaderToBed {
    is => ['Command'],
    has_input => [
        bam_file => {
            is => 'Text',
            doc => 'The BAM file to create BED file from.',
        },
        bed_file => {
            is => 'Text',
            doc => 'The output BED file of regions based on the provided BAM header.',
            is_output => 1,
        },
    ],
};

sub help_brief {
    "A BAM header can be converted to BED coordinates for each sequence in the sequence dictionary.  (THIS CAN BE REPLACED BY USING A BAM AS THE ROI FILE!)",
}
sub help_detail {
return <<EOS
A BAM file contains a sequence dictionary with the name and length of each sequence
in the reference aligned too.  This command generates a BED format file with
each reference sequence as an ROI.  This is no longer necessary since ref-cov
will now accept BAM files as an roi-file-path.  Using the BAM file as an
roi-file-path or generating a BED file to be used as ROI allows for coverage
evaluations across entire chromosomes.
EOS
}
sub execute {
    my $self = shift;

    my $bam = Genome::Model::Tools::RefCov::Bam->create(
        bam_file => $self->bam_file,
    );
    my @headers = qw/chr start end name/;
    my $bed_writer = Genome::Utility::IO::SeparatedValueWriter->create(
        output => $self->bed_file,
        separator => "\t",
        headers => \@headers,
        print_headers => 0,
    );
    my $chr_to_length = $bam->chr_to_length_hash_ref;
    # TODO: Replace with a sorting algorithm that matches samtools/picard
    for my $chr (sort {$chr_to_length->{$b} <=> $chr_to_length->{$a}} keys %{$chr_to_length}) {
        my %data;
        $data{chr} = $chr;
        $data{start} = 0;
        $data{end} = $chr_to_length->{$chr};
        $data{name} = $chr;
        $bed_writer->write_one(\%data);
    }
    return 1;
}
