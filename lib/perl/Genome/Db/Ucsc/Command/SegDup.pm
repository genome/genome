package Genome::Db::Ucsc::Command::SegDup;

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
    "otherChrom",
    "otherStart",
    "otherEnd",
    "otherSize",
    "uid",
    "posBasesHit",
    "testResult",
    "verdict",
    "chits",
    "ccov",
    "alignfile",
    "alignL",
    "indelN",
    "indelS",
    "alignB",
    "matchB",
    "mismatchB",
    "transitionsB",
    "transversionsB",
    "fracMatch",
    "fracMatchIndel",
    "jcK",
    "k2K",
);

class Genome::Db::Ucsc::Command::SegDup {
    is => 'Genome::Db::Ucsc::Command::Base',
    doc => "Fetches gap data for a reference name",
    has => [
        sort_position => {
            default_value => '2,4',
        },
    ],
};

sub headings {
    my $self = shift;
    return \@_HEADINGS;
}

sub table_names {
    my $self = shift;
    return ['genomicSuperDups'];
}

1;
