package Genome::VariantReporting::Vep::SplitAlternateAlleles;

use strict;
use warnings FATAL => 'all';
use Genome;
use Genome::File::Vcf::Reader;

class Genome::VariantReporting::Vep::SplitAlternateAlleles {
    is => 'Command::V2',
    has_input => [
        input_file => {
            is => 'Path',
        },
    ],
    has_output => [
        output_file => {
            is => 'Path',
        },
    ],
};

my @COLUMN_HEADERS = qw(
    CHROM
    POS
    ID
    REF
    ALT
    QUAL
    FILTER
    INFO
);

sub execute {
    my $self = shift;

    my $reader = Genome::File::Vcf::Reader->new($self->input_file);
    my $outfile = Genome::Sys->open_file_for_writing($self->output_file);

    print $outfile get_header($reader->header);

    while (my $entry = $reader->next()) {
        my @alternate_alleles = @{$entry->{alternate_alleles}};
        for my $alt (@alternate_alleles) {
            print $outfile get_line($entry, $alt);
        }
    }
    $outfile->close();

    return 1;
}

sub get_line {
    my ($entry, $alt) = @_;

    return join("\t",
        $entry->{chrom},
        $entry->{position},
        join(",", @{$entry->{identifiers}}) || '.',
        $entry->{reference_allele} || '.',
        $alt,
        $entry->{quality} || '.',
        join(";", @{$entry->{_filter}}) || '.',
        '.',
    ) . "\n";
}

sub get_header {
    my $header = shift;

    my @unparsed = map {"##$_"} @{$header->_unparsed_lines};

    my @lines;
    push @lines, sprintf("##fileformat=%s", $header->fileformat) if defined $header->fileformat;
    push @lines, grep {defined} (
        @unparsed,
        #$header->_metainfo_lines,
        #$header->_info_lines,
        $header->_format_lines,
        $header->_filter_lines,
        #$header->_header_line
        "#" . join("\t", @COLUMN_HEADERS),
        );
    return join("\n", @lines) . "\n";
}

1;
