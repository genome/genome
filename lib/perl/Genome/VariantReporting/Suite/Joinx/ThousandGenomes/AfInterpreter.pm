package Genome::VariantReporting::Suite::Joinx::ThousandGenomes::AfInterpreter;

use strict;
use warnings;
use Genome;

class Genome::VariantReporting::Suite::Joinx::ThousandGenomes::AfInterpreter {
    is => ['Genome::VariantReporting::Framework::Component::Interpreter'],
    doc => 'Output the allele frequency based on 1000Genomes',
};

sub name {
    return '1kg'
}

sub requires_annotations {
    return qw/1kg/;
}

sub field_descriptions {
    return (
        '1kg-af' => 'Allele frequency based on thousand genomes project'
    );
}

sub _interpret_entry {
    my $self = shift;
    my $entry = shift;
    my $passed_alt_alleles = shift;

    my %return_values;
    my $af = $entry->info("AF");
    for my $variant_allele (@$passed_alt_alleles) {
        $return_values{$variant_allele}->{'1kg-af'} = $af;
    }
    return %return_values;
}
1;

