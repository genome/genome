package Genome::Annotation::FTKeepFilter;

use strict;
use warnings FATAL => 'all';
use Genome;
use List::Util qw/first/;

class Genome::Annotation::FTKeepFilter {
    is => 'Genome::Annotation::FilterBase',
    has => {
        keep_filter_values => {
            is => 'Text',
            is_many => 1,
            default_value => ['PASS'],
        },
        sample_index => {
            is => 'Integer',
        },
    },
};

sub name {
    return 'ft-keep';
}

sub process_entry {
    my ($self, $entry) = @_;

    my $ft_string = $entry->sample_field($self->sample_index, 'FT');
    return $self->return_values($entry, 0) unless defined($ft_string);

    my @ft_values = split(';', $ft_string);
    for my $ft_value (@ft_values) {
        if (first { $_ eq $ft_value } $self->keep_filter_values) {
            return $self->return_values($entry, 1);
        }
    }

    return $self->return_values($entry, 0);
}

sub return_values {
    my $self = shift;
    my $entry = shift;
    my $pass = shift;

    my %return_values;
    for my $alt_allele (@{$entry->{alternate_alleles}}) {
        $return_values{$alt_allele} = $pass;
    }
    return %return_values;
}

1;
