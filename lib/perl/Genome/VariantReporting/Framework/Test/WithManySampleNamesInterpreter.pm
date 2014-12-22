package Genome::VariantReporting::Framework::Test::WithManySampleNamesInterpreter;

use strict;
use warnings FATAL => 'all';
use Genome;
use Data::Dump qw(pp);

class Genome::VariantReporting::Framework::Test::WithManySampleNamesInterpreter {
    is => ['Genome::VariantReporting::Framework::Component::Interpreter', 'Genome::VariantReporting::Framework::Component::WithManySampleNames'],
};

sub name {
    return '__with_many_sample_names__';
}

sub requires_annotations {
    return qw(__test__);
}

sub available_fields {
    my $self = shift;
    my %field_descriptions = $self->field_descriptions;
    return $self->create_sample_specific_field_names([keys %field_descriptions]);
}

sub field_descriptions {
    return (
        info => "info description",
    );
}

sub interpret_entry {
    my $self = shift;
    my $entry = shift;
    my $passed_alt_alleles = shift;

    my %return_values;
    for my $variant_allele (@$passed_alt_alleles) {
        for my $field ($self->available_fields) {
            $return_values{$variant_allele}->{$field} = pp($entry->info_for_allele($variant_allele));
        }
    }

    return %return_values;
}

1;
