package Genome::VariantReporting::Generic::MaxInfoThresholdFilter;

use strict;
use warnings;
use Genome;

class Genome::VariantReporting::Generic::MaxInfoThresholdFilter {
    is  => ['Genome::VariantReporting::Framework::Component::Filter'],
    has => [
        info_tag => {
            is  => 'String',
            doc => 'custom tag name in the info field to show threshold, like SEGDUP=82',
        },
        threshold => {
            is  => 'Number',
            doc => 'The maximum thresold to pass',
        }
    ],
    doc => 'Filter out variants that exceed the specified threshold in the specified INFO field',
};

sub name {
    return 'max-info-threshold';
}

sub requires_annotations {
    return ();
}

sub filter_entry {
    my ($self, $entry) = @_;
    my $threshold = $self->threshold;

    unless (defined $threshold and $threshold =~ /^\d+\.?\d*$/) {
        die "the input threshold: $threshold is not a number\n";
    }

    unless (exists $entry->{header}->info_types->{$self->info_tag}) {
        die sprintf("VCF header does not contain info tag (%s)", $self->info_tag);
    }

    my %return_values;

    for my $alt_allele ( @{$entry->{alternate_alleles}} ) {
        my $number = $entry->info($self->info_tag);
        $return_values{$alt_allele} = $number > $threshold ? 0 : 1;
    }

    return %return_values;
}

sub vcf_id {
    my $self = shift;
    return $self->info_tag , '_' . $self->threshold;
}

sub vcf_description {
    my $self = shift;
    return 'Filter out variants with info field ' . $self->info_tag . ' being greater than ' . $self->threshold;
}

1;
