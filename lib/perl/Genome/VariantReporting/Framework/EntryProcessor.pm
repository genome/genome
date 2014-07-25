package Genome::VariantReporting::Framework::EntryProcessor;

use strict;
use warnings;
use Genome;

class Genome::VariantReporting::Framework::EntryProcessor {
    has => [
        reporter => {
            is => 'Genome::VariantReporting::Framework::Component::Reporter',
        },
        filters => {
            is => 'ARRAY',
        },
        interpreters => {
            is => 'ARRAY',
        },
    ],
};

sub process_entry {
    my $self = shift;
    my $entry = shift;

    my $interpretations = $self->interpretations($entry);
    return unless keys %$interpretations;

    $self->reporter->report($interpretations);
}

sub interpretations {
    my $self = shift;
    my $entry = shift;

    my %interpretations;
    my $passed_alleles = $self->passed_alleles($entry);
    return {} unless scalar(@$passed_alleles);

    for my $interpreter (@{$self->interpreters}) {
        $interpretations{$interpreter->name} = {$interpreter->interpret_entry($entry, $passed_alleles)};
    }

    return \%interpretations;
}

sub passed_alleles {
    my $self = shift;
    my $entry = shift;

    my $filter_results = initialize_filters($entry);
    for my $filter (@{$self->filters}) {
        combine($filter_results, {$filter->filter_entry($entry)});
        last if(all_zeros($filter_results));
    }

    return [grep {$filter_results->{$_} == 1} keys %$filter_results];
}

sub initialize_filters {
    my $entry = shift;
    my %filter_values;
    for my $allele (@{$entry->{alternate_alleles}}) {
        $filter_values{$allele} = 1;
    }
    return \%filter_values;
}

sub combine {
    my $accumulator = shift;
    my $new_result = shift;
    for my $allele (keys %$accumulator) {
        $accumulator->{$allele} = $accumulator->{$allele} & $new_result->{$allele};
    }
    return $accumulator;
}

sub all_zeros {
    my $filter_results = shift;

    for my $key (keys %$filter_results) {
        return 0 if $filter_results->{$key};
    }

    return 1;
}

1;

