package Genome::Db::Ucsc::Command::GapList;

use strict;
use warnings;

use Genome;

our @_HEADINGS = (
    "chrom",
    "chromStart",
    "chromEnd",
    "ix",
    "n",
    "size",
    "type",
    "bridge",
    "bin",
);

our %_DB_TABLES_NAMES = (
    'hg38' => ['gap'],
    'hg19' => ['gap'],
    'hg18' => [
        'chr1_gap',
        'chr2_gap',
        'chr2_gap',
        'chr3_gap',
        'chr4_gap',
        'chr5_gap',
        'chr6_gap',
        'chr7_gap',
        'chr8_gap',
        'chr9_gap',
        'chr10_gap',
        'chr11_gap',
        'chr12_gap',
        'chr13_gap',
        'chr14_gap',
        'chr15_gap',
        'chr16_gap',
        'chr17_gap',
        'chr18_gap',
        'chr19_gap',
        'chr20_gap',
        'chr21_gap',
        'chr22_gap',
        'chrX_gap',
        'chrY_gap',
    ],
);


class Genome::Db::Ucsc::Command::GapList {
    is => 'Genome::Db::Ucsc::Command::Base',
    doc => "Fetches gap data for a reference name",
    has => [
        sort_position => {
            default_value => '1,3',
        },
    ],
};

sub headings {
    my $self = shift;
    return \@_HEADINGS;
}

sub table_names {
    my $self = shift;
    return $_DB_TABLES_NAMES{$self->reference_name};
}


1;
