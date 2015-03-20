package Genome::VariantReporting::Generic::VariantCallersInterpreter;

use strict;
use warnings;
use Genome;
use List::MoreUtils qw/uniq/;

class Genome::VariantReporting::Generic::VariantCallersInterpreter {
    is => ['Genome::VariantReporting::Framework::Component::Interpreter', 'Genome::VariantReporting::Framework::Component::WithSampleName'],
    doc => 'Output a list of variant callers that called this position for the specified sample, as well as the number of callers',
    has => [
        valid_callers => {
            is => 'String',
            is_many => 1,
            doc => 'List of variant callers to include in determination for filtering',
        },
    ],
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
    my %callers;
    for my $alt_allele (@$passed_alt_alleles) {
        $return_values{$alt_allele} = { variant_callers => [] };
        $callers{$alt_allele} = [];
    }

    for my $caller_name ($self->valid_callers) {
        my $sample_name = $self->sample_name_with_suffix($caller_name);
        my $sample_index = eval{ $entry->{header}->index_for_sample_name($sample_name) };
        my $error = $@;
        if (defined($error) &&
            $error =~ /^\QSample name $sample_name not found in header\E/) {
            next;
        }
        my @sample_alt_alleles = $entry->alt_bases_for_sample($sample_index);
        for my $sample_alt_allele (uniq @sample_alt_alleles) {
            push(@{$callers{$sample_alt_allele}}, $caller_name);
        }
    }

    for my $alt_allele (keys %return_values) {
        $return_values{$alt_allele} = {
            variant_callers => $callers{$alt_allele},
            variant_caller_count => scalar(@{$callers{$alt_allele}}),
        };
    }

    return %return_values;
}
1;

