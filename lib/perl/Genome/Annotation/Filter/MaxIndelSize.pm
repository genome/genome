package Genome::Annotation::Filter::MaxIndelSize;

use strict;
use warnings;
use Genome;

class Genome::Annotation::Filter::MaxIndelSize {
    is => 'Genome::Annotation::Filter::Base',
    has => [
        size => {
            is => "Number",
            doc => "The mazimum size of an INDEL to pass",
        },
    ],
};

sub __errors__ {
    my $self = shift;

    my @errors = $self->SUPER::__errors__;
    return @errors if @errors;

    my $size = $self->size;
    if ( $size !~ /^\d+$/ ) {
        push @errors, UR::Object::Tag->create(
            type => 'error',
            properties => [qw/ size /],
            desc => "Value given ($size) is not a whole number!",
        );
    }

    return @errors;
}

sub name {
    return 'indel-size';
}

sub filter_entry {
    my $self = shift;
    my $entry = shift;

    my %return_values;
    for my $alt_allele ( @{$entry->{alternate_alleles}} ) {
        my $indel_legnth = abs( length($entry->{reference_allele}) - length($alt_allele) );
        $return_values{$alt_allele} = ( $indel_legnth >= $self->size ) ? 1 : 0;
    }

    return %return_values;
}

1;

