package Genome::VariantReporting::Expert::BamReadcount::MinCoverageObservedInterpreter;

use strict;
use warnings;
use Genome;
use List::Util qw/ min /;
use Scalar::Util qw( looks_like_number );

class Genome::VariantReporting::Expert::BamReadcount::MinCoverageObservedInterpreter {
    is => ['Genome::VariantReporting::Framework::Component::Interpreter', 'Genome::VariantReporting::Framework::Component::WithManySampleNames'],
    has => [],
};

sub name {
    return 'min-coverage-observed';
}

sub requires_annotations {
    return ('bam-readcount');
}

sub field_descriptions {
    return (
        min_coverage_observed => 'Minimum coverage (ref_count+var_count) between all the samples at this position'
    );
}

sub _interpret_entry {
    my $self = shift;
    my $entry = shift;
    my $passed_alt_alleles = shift;

    my %coverages;
    for my $sample_name ($self->sample_names) {
        my $interpreter = Genome::VariantReporting::Expert::BamReadcount::VafInterpreter->create(sample_name => $sample_name);
        my %result = $interpreter->interpret_entry($entry, $passed_alt_alleles);
        for my $alt_allele (@$passed_alt_alleles) {
            $coverages{$alt_allele}->{coverage}->{$sample_name} = $result{$alt_allele}->{ref_count} + $result{$alt_allele}->{var_count};
        }
    }

    my %return_values;
    for my $alt_allele (@$passed_alt_alleles) {
        my @coverages = grep { looks_like_number($_) } values %{$coverages{$alt_allele}->{coverage}};
        if (@coverages) {
            $return_values{$alt_allele} = { min_coverage_observed => min(@coverages) };
        } else {
            $return_values{$alt_allele} = { min_coverage_observed => $self->interpretation_null_character };
        }
    }
    return %return_values;
}

1;
