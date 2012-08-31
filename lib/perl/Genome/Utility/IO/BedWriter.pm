package Genome::Utility::IO::BedWriter;

use strict;
use warnings;

use Genome;

class Genome::Utility::IO::BedWriter {
    is => ['Genome::Utility::IO::SeparatedValueWriter'],
    has => [
        separator => {
            default_value => "\t",
        },
        print_headers => {
            default_value => 0,
        },
        ignore_extra_columns => {
            default_value => 1,
        },
    ],
};


sub create {
    my $class = shift;
    my %params = @_;

    my $headers = delete $params{headers};
    unless ($headers) {
        $headers = Genome::Utility::IO::BedReader->default_headers;
    }
    $params{headers} = $headers;
    my $self = $class->SUPER::create(%params)
        or return;
    return $self;
}

1;
