package Genome::VariantReporting::Framework::Component::Interpreter::Fpkm;

use strict;
use warnings;
use Genome;

class Genome::VariantReporting::Framework::Component::Interpreter::Fpkm {
    is => ['Genome::VariantReporting::Framework::Component::WithSampleName'],
};

sub fpkm_for_genotype_allele {
    my ($self, $entry) = @_;

    my $sample_index = $self->sample_index($entry->{header});
    my @genotype_alleles = $entry->bases_for_sample($sample_index);

    my $fpkm = $entry->sample_field($self->sample_index($entry->{header}), "FPKM");
    my @fpkms;
    if ($fpkm) {
        @fpkms = split(',', $fpkm);
    }
    else {
        @fpkms = ('.') x @genotype_alleles;
    }

    if (scalar(@fpkms) ne scalar(@genotype_alleles)) {
        die $self->error_message("There should be the same number of FPKM values (%s) as there are genotype alleles (%s).", join(',', @fpkms), join(',', @genotype_alleles));
    }

    my %fpkm_for_genotype_allele;
    for my $i (0..$#fpkms) {
        $fpkm_for_genotype_allele{$genotype_alleles[$i]} = $fpkms[$i];
    }

    return %fpkm_for_genotype_allele;
}

1;
