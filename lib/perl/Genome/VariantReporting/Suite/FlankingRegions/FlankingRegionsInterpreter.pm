package Genome::VariantReporting::Suite::FlankingRegions::FlankingRegionsInterpreter;

use strict;
use warnings;
use Genome;

class Genome::VariantReporting::Suite::FlankingRegions::FlankingRegionsInterpreter {
    is => ['Genome::VariantReporting::Framework::Component::Interpreter'],
    doc => 'Output the flanking sequence for the reference (FSRF info field) and the alternate allele (FSAF info field)',
};

sub name {
    return 'flanking-regions';
}

sub requires_annotations {
    return ('flanking-regions');
}

sub field_descriptions {
    return (
        reference_fasta => 'Reference allele with flanking region',
        alt_fasta => 'Alternate allele with flanking region',
    );
}

sub interpret_entry {
    my $self = shift;
    my $entry = shift;
    my $passed_alt_alleles = shift;

    my %return_values;
    for my $alt_allele (@$passed_alt_alleles) {
        $return_values{$alt_allele} = {
            alt_fasta => $entry->info_for_allele($alt_allele, 'FSAF'),
            reference_fasta => $entry->info->{FSRF},
        }
    }
    return %return_values;
}

1;

