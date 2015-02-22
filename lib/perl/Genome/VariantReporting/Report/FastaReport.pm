package Genome::VariantReporting::Report::FastaReport;

use strict;
use warnings;
use Genome;

class Genome::VariantReporting::Report::FastaReport {
    is => 'Genome::VariantReporting::Framework::Component::Report::SingleFile',
    has_transient_optional => [
        _output_fh => {
            is_structural => 1,
        },
    ],
    doc => 'Prints out a fasta file with both reference and alternate alleles surrounded by a short flanking sequence',
};

sub name {
    return 'fasta';
}

sub required_interpreters {
    return qw(position flanking-regions vaf vep);
}

sub report {
    my $self = shift;
    my $interpretations = shift;
    my $long_allele_counter = 0;
    for my $allele (keys %{$interpretations->{($self->requires_interpreters)[0]}}) {
        my $header_allele = $allele;
        if (length($allele) > 50) {
            $header_allele = "long_allele$long_allele_counter";
            $long_allele_counter++;
        }
        my $fasta_header = sprintf(">%s:%s:%s:%s:%s:%s",
            $interpretations->{position}->{$allele}->{chromosome_name},
            $interpretations->{position}->{$allele}->{start},
            $header_allele,
            $interpretations->{'vaf'}->{$allele}->{vaf},
            $interpretations->{'vep'}->{$allele}->{default_gene_name},
            $interpretations->{'vep'}->{$allele}->{ensembl_gene_id});
        $self->_output_fh->print($fasta_header."-mut\n");
        $self->_output_fh->print(
            $interpretations->{'flanking-regions'}->{$allele}->{alt_fasta}."\n");
        $self->_output_fh->print($fasta_header."-wt\n");
        $self->_output_fh->print(
            $interpretations->{'flanking-regions'}->{$allele}->{reference_fasta}."\n");
    }
}


1;
