package Genome::VariantReporting::Generic::VcfEntryInterpreter;

use strict;
use warnings FATAL => 'all';
use Genome;
use Set::Scalar;

class Genome::VariantReporting::Generic::VcfEntryInterpreter {
    is => ['Genome::VariantReporting::Framework::Component::Interpreter'],
};

sub name {
    return 'vcf-entry';
}

sub requires_annotations {
    ();
}

sub field_descriptions {
    return (
        vcf_entry => 'Complete vcf entry object',
    );
}

sub _interpret_entry {
    my ($self, $entry, $passed_alt_alleles) = @_;

    my %return_values;
    for my $alt_allele (@$passed_alt_alleles) {
        $return_values{$alt_allele} =  {
            vcf_entry => $entry
        };
    }
    return %return_values;
}


1;
