package Genome::Utility::IO::GenePredReader;

use strict;
use warnings;

use Genome;

my @DEFAULT_HEADERS = qw/name chrom strand txStart txEnd cdsStart cdsEnd exonCount exonStarts exonEnds id name2 cdsStartStat cdsEndStat exonFrames/;

class Genome::Utility::IO::GenePredReader {
    is => ['Genome::Utility::IO::SeparatedValueReader'],
    has => [
        separator => {
            default_value => "\t",
        },
    ],
};

sub default_headers {
    return \@DEFAULT_HEADERS;
}

sub create {
    my $class = shift;
    my %params = @_;

    my $headers = delete $params{headers};
    unless ($headers) {
        $headers = $class->default_headers;
    }
    $params{headers} = $headers;
    my $self = $class->SUPER::create(%params)
        or return;
    return $self;
}


1;
