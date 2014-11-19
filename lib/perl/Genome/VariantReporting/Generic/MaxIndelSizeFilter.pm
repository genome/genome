package Genome::VariantReporting::Generic::MaxIndelSizeFilter;

use strict;
use warnings;
use Genome;

class Genome::VariantReporting::Generic::MaxIndelSizeFilter {
    is => 'Genome::VariantReporting::Framework::Component::Filter',
    has => [
        size => {
            is => "Number",
            doc => "The maximum size of an INDEL to pass",
        },
    ],
    doc => 'Filter out indels that exceed the specified size',
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
    return 'max-indel-size';
}

sub requires_annotations {
    return qw/ /;
}

sub filter_entry {
    my $self = shift;
    my $entry = shift;

    my %return_values;
    for my $alt_allele ( @{$entry->{alternate_alleles}} ) {
        my $indel_length = abs( length($entry->{reference_allele}) - length($alt_allele) );
        $return_values{$alt_allele} = ( $indel_length >= $self->size ) ? 0 : 1;
    }

    return %return_values;
}

sub vcf_id {
    my $self = shift;
    return 'MAXINDEL' . $self->size;
}

sub vcf_description {
    my $self = shift;

    return 'The size of an indel is less than or equal to ' . $self->size;
}

1;
