package Genome::VariantReporting::Generic::FTKeepFilter;

use strict;
use warnings FATAL => 'all';
use Genome;
use List::Util qw/first/;

class Genome::VariantReporting::Generic::FTKeepFilter {
    is => ['Genome::VariantReporting::Framework::Component::Filter', 'Genome::VariantReporting::Framework::Component::WithSampleName'],
    has => {
        keep_filter_values => {
            is => 'Text',
            is_many => 1,
            default_value => ['PASS'],
            doc => 'A list of FT values',
        },
    },
    doc => q{Filter out variants that don't contain the specified FT values},
};

sub name {
    return 'ft-keep';
}

sub requires_annotations {
    ();
}

sub filter_entry {
    my ($self, $entry) = @_;

    my $ft_string = $entry->sample_field($self->sample_index($entry->{header}), 'FT');
    return $self->pass_return_values($entry) unless defined($ft_string);

    my @ft_values = split(';', $ft_string);
    for my $ft_value (@ft_values) {
        if (first { $_ eq $ft_value } $self->keep_filter_values) {
            return $self->pass_return_values($entry);
        }
    }

    return fail_return_values($entry);
}

sub pass_return_values {
    my ($self, $entry) = @_;

    my %return_values = fail_return_values($entry);
    my @sample_alt_alleles = $entry->alt_bases_for_sample($self->sample_index($entry->{header}));
    for my $alt_allele (@sample_alt_alleles) {
        $return_values{$alt_allele} = 1;
    }
    return %return_values;
}

sub fail_return_values {
    my $entry = shift;

    return map { $_ => 0 } @{$entry->{alternate_alleles}};
}

sub vcf_id {
    my $self = shift;
    return 'FT' . ($self->keep_filter_values)[0];
}

sub vcf_description {
    my $self = shift;
    return sprintf(
        'FT value for sample %s is one of [%s]',
        $self->sample_name,
        join(",", $self->keep_filter_values)
    );
}

1;
