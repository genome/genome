package Genome::VariantReporting::Generic::GenotypeInterpreter;

use strict;
use warnings;

use Genome;
use Set::Scalar;

class Genome::VariantReporting::Generic::GenotypeInterpreter {
      is => ['Genome::VariantReporting::Framework::Component::Interpreter', 'Genome::VariantReporting::Framework::Component::WithSampleName'],
};

sub name {
    return 'genotype';
}

sub requires_annotations {
    return qw/ /;
}

sub field_descriptions {
    my $self = shift;
    return (
        genotype => sprintf('Genotype for sample %s', $self->sample_name),
    );
}

sub _interpret_entry {
    my $self = shift;
    my $entry = shift;
    my $passed_alt_alleles = shift;

    my $sample_alt_alleles = Set::Scalar->new($entry->alt_bases_for_sample($self->sample_index($entry->{header})));
    my $gt = $entry->genotype_for_sample($self->sample_index($entry->{header}));
    unless ($sample_alt_alleles and $gt) {
        return $self->null_interpretation($passed_alt_alleles);
    }

    my %return_values;
    for my $alt_allele ( @$passed_alt_alleles ) {
        unless ($sample_alt_alleles->contains($alt_allele)) {
            $return_values{$alt_allele} =  { genotype => 'not called' };
            next;
        }

        if ( $gt->is_missing ) {
            # fail
            $return_values{$alt_allele} =  { genotype => 'missing' };
        }
        elsif ( $gt->is_homozygous ) {
            $return_values{$alt_allele} = { genotype => 'homozygous' };
        }
        else { # $gt->is_heterozygous
            $return_values{$alt_allele} = { genotype => 'heterozygous' };
        }
    }

    return %return_values;
}

1;

