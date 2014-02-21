package Genome::Model::GenotypeMicroarray::GenotypeFile::DefaultHeader;

use strict;
use warnings;

use Genome;

use Data::Dumper 'Dumper';

class Genome::Model::GenotypeMicroarray::GenotypeFile::DefaultHeader { 
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
    {
        name => 'alleles', # origninal_genotype
        id => 'OG',
        header => ',Number=1,Type=Integer,Description="Original genotype calls">', 
    },
];
sub supported_info_fields {
    return $supported_info_fields;
}

my $header_lines;
sub header_lines {

    return $header_lines if $header_lines;

    $header_lines = [
        '##fileformat=VCFv4.1',
    ];
    for my $field ( @$supported_info_fields ) {
        push @$header_lines, '##INFO=<ID='.$field->{id}.$field->{header};
    }
    push @$header_lines, '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">';

    return $header_lines;
}

our $header;
sub header {
    return $header if $header;
    $header = Genome::File::Vcf::Header->create(lines => header_lines());
    return $header;
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

1;

