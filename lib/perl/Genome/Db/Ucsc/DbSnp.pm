package Genome::Db::Ucsc::DbSnp;

use strict;
use warnings;

use Genome;

our @_HEADINGS = (
    "bin",
    "chrom",
    "chromStart",
    "chromEnd",
    "name",
    "score",
    "strand",
    "refNCBI",
    "observed",
    "molType",
    "class",
    "valid",
    "avHet",
    "avHetSE",
    "func",
    "locType",
    "weight",
);

our %_DB_TABLES_NAMES = (
    "137" => {'hg19' => ['snp137'],
    'hg18' => [],}
);


class Genome::Db::Ucsc::DbSnp {
    is => 'Genome::Db::Ucsc::Base',

    doc => "Fetches gap data for a reference name",
    has_input => [
        dbsnp_version => {
            is => 'Text',
        },
    ],
};

sub headings {
    my $self = shift;
    return \@_HEADINGS;
}

sub table_names {
    my $self = shift;
    return $_DB_TABLES_NAMES{$self->dbsnp_version}->{$self->reference_name};
}

sub filter {
    my $self = shift;
    return "class='named'";
}

1;
