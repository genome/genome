package Genome::VariantReporting::BamReadcount::MinCoverageObservedInterpreter;

use strict;
use warnings;
use Genome;
use List::Util qw/ min /;

class Genome::VariantReporting::BamReadcount::MinCoverageObservedInterpreter {
    is => ['Genome::VariantReporting::Framework::Component::Interpreter', 'Genome::VariantReporting::Framework::Component::WithManySampleNames'],
    has => [],
};

sub name {
    return 'min-coverage-observed';
}

sub requires_annotations {
    return ('bam-readcount');
}


sub available_fields {
    return qw/min_coverage_observed/;
}

sub interpret_entry {
    my $self = shift;
    my $entry = shift;
    my $passed_alt_alleles = shift;

    my %coverages;
    for my $sample_name ($self->sample_names) {
        my $interpreter = Genome::VariantReporting::BamReadcount::VafInterpreter->create(sample_name => $sample_name);
        my %result = $interpreter->interpret_entry($entry, $passed_alt_alleles);
        for my $alt_allele (@$passed_alt_alleles) {
            $coverages{$alt_allele}->{ref_count}->{$sample_name} = $result{$alt_allele}->{ref_count};
            $coverages{$alt_allele}->{var_count}->{$sample_name} = $result{$alt_allele}->{var_count};
        }
    }

    my %return_values;
    for my $alt_allele (@$passed_alt_alleles) {
        my @ref_counts = values %{$coverages{$alt_allele}->{ref_count}};
        my @var_counts = values %{$coverages{$alt_allele}->{var_count}};
        $return_values{$alt_allele} = { min_coverage_observed => min(@ref_counts, @var_counts) };
    }
    return %return_values;
}

1;
