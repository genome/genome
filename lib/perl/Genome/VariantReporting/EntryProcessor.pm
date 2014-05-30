package Genome::VariantReporting::EntryProcessor;

use strict;
use warnings;
use Genome;
use Memoize qw();

class Genome::VariantReporting::EntryProcessor {
    has => [
        reporter_plan => {
            is => 'Genome::VariantReporting::Plan::ReporterPlan',
        },
        translations => {
            is => 'HASH'
        },
    ],
};

sub process_entry {
    my $self = shift;
    my $entry = shift;
    my $reporter_plan = shift;

    my $interpretations = $self->interpretations($entry);
    $self->report->report($interpretations);
}

sub report {
    my $self = shift;
    return $self->reporter_plan->object;
}

sub interpretations {
    my $self = shift;
    my $entry = shift;

    my %interpretations;
    my $passed_alleles = $self->passed_alleles($entry);
    for my $interpreter ($self->interpreters) {
        $interpretations{$interpreter->name} = {$interpreter->interpret_entry($entry, $passed_alleles)};
    }

    return \%interpretations;
}

sub interpreters {
    my $self = shift;
    my @interpreters;
    for my $plan ($self->reporter_plan->interpreter_plans) {
        push @interpreters, $plan->object($self->_translated_params($plan));
    }
    return @interpreters;
}
Memoize::memoize('interpreters');

sub passed_alleles {
    my $self = shift;
    my $entry = shift;

    my $filter_results = initialize_filters($entry);
    for my $filter ($self->filters) {
        combine($filter_results, {$filter->filter_entry($entry)});
        last if(all_zeros($filter_results));
    }

    return [grep {$filter_results->{$_} == 1} keys %$filter_results];
}

sub filters {
    my $self = shift;
    my @filters;
    for my $plan ($self->reporter_plan->filter_plans) {
        push @filters, $plan->object($self->_translated_params($plan));
    }
    return @filters;
}
Memoize::memoize('filters');

sub _translated_params {
    my ($self, $plan) = @_;
    my %params;
    for my $input_name ($plan->get_class->translated_inputs) {
        my $translated_value = $self->translations->{$plan->params->{$input_name}};
        if (defined($translated_value)) {
            $params{$input_name} = $translated_value;
        } else {
            die $self->error_message("Cannot translate input (%s) for filter (%s)",
                $input_name, $plan->name);
        }
    }
    return %params;
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

