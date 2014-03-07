package Genome::Model::GenotypeMicroarray::GenotypeFile::VcfHelper;

use strict;
use warnings;

use Genome;

class Genome::Model::GenotypeMicroarray::GenotypeFile::VcfHelper { 
    is => 'UR::Singleton',
};

our $supported_info_fields = [
    # SORTED BY ID!
    {
        name => 'allele_frequency',
        id => 'AF',
        header => ',Number=1,Type=Float,Description="Allele frequency for each ALT allele in the same order as listed">',
    },
    {
        name => 'cnv_confidence',
        id => 'CC',
        header => ',Number=1,Type=Float,Description="CNV Confidence">',
    },
    {
        name => 'cnv_value',
        id => 'CV',
        header => ',Number=1,Type=Float,Description="CNV Value">',
    },
    {
        name => 'log_r_ratio',
        id => 'LR',
        header => ',Number=1,Type=Float,Description="Log R Ratio">',
    },
    {
        name => 'gc_score',
        id => 'GC',
        header => ',Number=1,Type=Float,Description="GC Score">',
    },
];
sub supported_info_fields {
    return $supported_info_fields;
}

our $info_order;
sub info_order { 
    return $info_order if $info_order;
    $info_order = [ map { $_->{id} } sort { $a->{id} cmp $b->{id} } @$supported_info_fields ];
    return $info_order;
}

our $info_names_and_ids;
sub info_names_and_ids { 
    return $info_names_and_ids if $info_names_and_ids;
    $info_names_and_ids = { map { $_->{name} => $_->{id} } @$supported_info_fields };
    return $info_names_and_ids;
}

our $formats = [qw/ GT OG /];
sub formats {
    return $formats
}

our $format_key_idx = { GT => 0, OG => 1, };
sub format_key_idx {
    return $format_key_idx;
}

sub header {
    my ($self, $sample_name) = @_;

    Carp::confess('No sample name given to get vcf header!') if not $sample_name;

    my @header_lines = (
        '##fileformat=VCFv4.1',
    );
    for my $field ( @$supported_info_fields ) {
        push @header_lines, '##INFO=<ID='.$field->{id}.$field->{header};
    }
    push @header_lines, '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">';
    push @header_lines, '##FORMAT=<ID=AL,Number=1,Type=String,Description="Alleles">';
    push @header_lines, '#'.join("\t", Genome::File::Vcf::Header->column_headers, $sample_name);

    my $header = Genome::File::Vcf::Header->create(lines => \@header_lines);
    return $header;
}

1;

