package Genome::VariantReporting::Generic::FTInterpreter;

use strict;
use warnings FATAL => 'all';
use Genome;
use Set::Scalar;

class Genome::VariantReporting::Generic::FTInterpreter {
    is => ['Genome::VariantReporting::Framework::Component::Interpreter', 'Genome::VariantReporting::Framework::Component::WithSampleName'],
    doc => 'Output the value of the FT field',
};

sub name {
    return 'ft';
}

sub requires_annotations {
    ();
}

sub field_descriptions {
    my $self = shift;
    return (
        ft_string => sprintf('FT sample field for sample %s', $self->sample_name),
    );
}

sub _interpret_entry {
    my ($self, $entry, $passed_alt_alleles) = @_;

    my $ft_string = $entry->sample_field($self->sample_index($entry->{header}), 'FT');

    my %return_values;

    my @sample_alt_alleles = $entry->alt_bases_for_sample($self->sample_index($entry->{header}));
    unless (@sample_alt_alleles) {
        return $self->null_interpretation($passed_alt_alleles);
    }

    my $sample_alt_alleles = Set::Scalar->new(@sample_alt_alleles);
    for my $alt_allele (@$passed_alt_alleles) {
        unless ($sample_alt_alleles->contains($alt_allele)) {
            $ft_string = "";
        }
        $return_values{$alt_allele} =  {
            ft_string => $ft_string
        };
    }
    return %return_values;
}


1;
