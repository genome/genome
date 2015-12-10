package Genome::VariantReporting::Report::EpitopeBindingPredictionReport;

use strict;
use warnings;
use Genome;
use List::AllUtils qw(all);

class Genome::VariantReporting::Report::EpitopeBindingPredictionReport {
    is => 'Genome::VariantReporting::Framework::Component::Report::SingleFile',
    doc => 'Output a fasta file of variant sequences for snps, frameshifts, and in frame indels for epitope binding prediction.'
};

sub name {
    return 'epitope-binding-prediction';
}

sub required_interpreters {
    return qw(epitope-variant-sequence);
}

sub allows_hard_filters {
    return 1;
}

sub file_name {
    return 'report.fasta';
}

sub report {
    my $self = shift;
    my $interpretations = shift;

    my %vcf_entry_interpretations = %{$interpretations->{'epitope-variant-sequence'}};

    for my $allele (keys %vcf_entry_interpretations) {
        my $variant_sequences = $vcf_entry_interpretations{$allele}->{variant_sequences};
        if ($variant_sequences) {
            while (my ($header, $sequence) = each %$variant_sequences) {
                $self->_output_fh->print("$header\n");
                $self->_output_fh->print("$sequence\n");
            }
        }
    }
}

1;
