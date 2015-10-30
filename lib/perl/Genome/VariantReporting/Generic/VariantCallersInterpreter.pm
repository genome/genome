package Genome::VariantReporting::Generic::VariantCallersInterpreter;

use strict;
use warnings;
use Genome;

class Genome::VariantReporting::Generic::VariantCallersInterpreter {
    is => ['Genome::VariantReporting::Framework::Component::Interpreter', 'Genome::VariantReporting::Generic::NCallersBase'],
    doc => 'Output a list of variant callers that called this position for the specified sample, as well as the number of callers',
};

sub name {
    return 'variant-callers';
}

sub requires_annotations {
    return ();
}

sub field_descriptions {
    my $self = shift;
    return (
        variant_callers => 'List of variant callers that called this position for sample ' . $self->sample_name,
        variant_caller_count => 'Number of callers that called this position for sample ' . $self->sample_name,
    );
}

sub _interpret_entry {
    my $self = shift;
    my $entry = shift;
    my $passed_alt_alleles = shift;

    my %return_values;
    for my $alt_allele (@$passed_alt_alleles) {
        $return_values{$alt_allele} = { variant_callers => [] };
    }

    my %callers = $self->get_callers($entry, $passed_alt_alleles);
    for my $alt_allele (keys %return_values) {
        $return_values{$alt_allele} = {
            variant_callers => $callers{$alt_allele},
            variant_caller_count => scalar(@{$callers{$alt_allele}}),
        };
    }

    return %return_values;
}
1;

