package Genome::Db::Ucsc::Command::RepeatMasker;

use strict;
use warnings;

use Genome;

our @_HEADINGS = (
    "bin",
    "genoName",
    "genoStart",
    "genoEnd",
    "repClass",
);

class Genome::Db::Ucsc::Command::RepeatMasker {
    is => 'Genome::Db::Ucsc::Command::Base',
    doc => "Fetches gap data for a reference name",
    has => [
        sort_position => {
            default_value => '2,4',
        },
    ]
};

sub headings {
    my $self = shift;
    return \@_HEADINGS;
}

sub table_names {
    my $self = shift;
    return ["rmsk"];
}

1;
