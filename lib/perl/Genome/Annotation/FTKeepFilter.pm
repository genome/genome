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
            is => 'Text',
        },
    },
};

sub name {
    return 'ft-keep';
}

sub requires_experts {
    return ();
}

sub process_entry {
    my ($self, $entry) = @_;

    my $ft_string = $entry->sample_field($self->sample_index, 'FT');
    return 0 unless defined($ft_string);

    my @ft_values = split(';', $ft_string);
    for my $ft_value (@ft_values) {
        if (first { $_ eq $ft_value } $self->keep_filter_values) {
            return 1;
        }
    }

    return 0;
}

1;
