package Genome::Utility::IO::GffReader;

use strict;
use warnings;

use Genome;

my @DEFAULT_HEADERS = qw/chr source type start end score strand frame attributes/;

class Genome::Utility::IO::GffReader {
    is => ['Genome::Utility::IO::SeparatedValueReader'],
    has => [
        ignore_lines_starting_with => {
            default_value => '##',
        },
        separator => {
            default_value => "\t",
        },
    ],
};

sub headers {
    return \@DEFAULT_HEADERS;
}

sub create {
    my $class = shift;
    my %params = @_;

    my $headers = delete $params{headers};
    unless ($headers) {
        $headers = $class->headers;
    }
    $params{headers} = $headers;
    my $self = $class->SUPER::create(%params)
        or return;
    return $self;
}

sub next_with_attributes {
    my $self = shift;

    my $data = $self->next;
    unless ($data) { return; }

    my $attributes = $data->{attributes};
    my @attributes = split(';',$attributes);
    my %attributes;
    for my $attribute (@attributes) {
        unless ($attribute =~ /^\s*(\S+)[\s=]*\s*\"?([^\s\"]+)\"?$/) {
            die(Data::Dumper::Dumper($data));
        }
        $data->{$1} = $2;
    }
    return $data;
}

sub next_with_attributes_hash_ref {
    my $self = shift;

    my $data = $self->next;
    unless ($data) { return; }

    my $attributes = $data->{attributes};
    my @attributes = split(';',$attributes);
    my %attributes;
    for my $attribute (@attributes) {
        unless ($attribute =~ /^\s*(\S+)[\s=]*\s*\"?([^\s\"]+)\"?$/) {
            die('Unknown attributes format: '. $attribute ."\n". Data::Dumper::Dumper($data));
        }
        $attributes{$1} = $2;
    }
    $data->{attributes_hash_ref} = \%attributes;
    return $data;
}

1;
