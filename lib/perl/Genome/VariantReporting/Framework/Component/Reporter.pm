package Genome::VariantReporting::Framework::Component::Reporter;

use strict;
use warnings FATAL => 'all';
use Genome;
use Carp qw(confess);
use List::MoreUtils qw(uniq);

class Genome::VariantReporting::Framework::Component::Reporter {
    is => [
        'Genome::VariantReporting::Framework::Component::WithTranslatedInputs',
    ],
    is_abstract => 1,
    has_transient_optional => [
        filters => {
            is => 'HASH',
            is_structural => 1,
            default => {},
        },
        interpreters => {
            is => 'HASH',
            is_structural => 1,
            default => {},
        },
    ],
};

sub name {
    my $self = shift;
    my $class = $self->class;
    confess "Abstract method 'name' must be defined in class '$class'";
}

sub required_interpreters {
    my $self = shift;
    my $class = $self->class;
    confess "Abstract method 'required_interpreters' must be defined in class '$class'";
}

sub allows_hard_filters {
    return 1;
}

sub initialize {
    # this gets called before interpretations are given to ->report method
    return;
}

sub report {
    my ($self, $interpretations) = shift;
    # do something with the interpretations to produce one or more reports
    return;
}

sub finalize {
    # this gets called after interpretations are given to ->report method
    return;
}

sub add_filter_object {
    my ($self, $filter) = @_;

    $self->filters->{$filter->name} = $filter;
}

sub add_filter_objects {
    my ($self, @filters) = @_;

    for my $filter (@filters) {
        $self->add_filter_object($filter);
    }
}

sub add_interpreter_object {
    my ($self, $interpreter) = @_;

    $self->interpreters->{$interpreter->name} = $interpreter;
}

sub add_interpreter_objects {
    my ($self, @interpreters) = @_;

    for my $interpreter (@interpreters) {
        $self->add_interpreter_object($interpreter);
    }
}

sub process_entry {
    my $self = shift;
    my $entry = shift;

    my $interpretations = $self->interpretations($entry);
    return unless keys %$interpretations;

    $self->report($interpretations);
}

sub interpretations {
    my $self = shift;
    my $entry = shift;

    my %interpretations;
    my $passed_alleles = $self->passed_alleles($entry);
    return {} unless scalar(@$passed_alleles);

    for my $interpreter (values %{$self->interpreters}) {
        $interpretations{$interpreter->name} = {
            $interpreter->interpret_entry($entry, $passed_alleles)
        };
    }

    return \%interpretations;
}

sub passed_alleles {
    my $self = shift;
    my $entry = shift;

    my $filter_results = initialize_filters($entry);
    for my $filter (values %{$self->filters}) {
        _boolean_combine($filter_results, {$filter->filter_entry($entry)});
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

sub _boolean_combine {
    my $accumulator = shift;
    my $new_result = shift;
    for my $allele (keys %$accumulator) {
        $accumulator->{$allele} = $accumulator->{$allele} & $new_result->{$allele};
    }
    return $accumulator;
}

sub all_zeros {
    my $filter_results = shift;

    for my $value (values %$filter_results) {
        return 0 if $value;
    }

    return 1;
}

sub vr_doc_sections {
    my $self = shift;

    my @sections = $self->SUPER::vr_doc_sections;

    if ($self->required_experts) {
        push @sections, {
            header => "REQUIRED EXPERTS",
            items => [$self->required_experts],
        };
    }
    if ($self->required_interpreters) {
        push @sections, {
            header => "REQUIRED INTERPRETERS",
            items => [$self->required_interpreters],
        };
    }
    return @sections;
}

sub required_experts {
    my $self = shift;

    my $factory = Genome::VariantReporting::Framework::Factory->create();
    my @experts;
    for my $interpreter ($self->required_interpreters) {
        my $class = $factory->get_class("interpreters", $interpreter);
        push @experts, $class->requires_annotations;
    }
    return uniq @experts;
}

1;
