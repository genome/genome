package Genome::Annotation::Interpreter::FpkmInterpreter;

use strict;
use warnings;
use Genome;

class Genome::Annotation::Interpreter::FpkmInterpreter {
    is => ['Genome::Annotation::Interpreter::Base', 'Genome::Annotation::WithSampleName'],
};

sub name {
    return 'fpkm';
}

sub requires_experts {
    ('fpkm');
}

sub available_fields {
    return qw/
        fpkm
    /;
}

sub interpret_entry {
    my $self = shift;
    my $entry = shift;
    my $passed_alt_alleles = shift;

    my %return_values;
    my %fpkm_for_genotype_allele = $self->fpkm_for_genotype_allele($entry);
    for my $variant_allele (@$passed_alt_alleles) {
        $return_values{$variant_allele} = {
            fpkm => $fpkm_for_genotype_allele{$variant_allele},
        };
    }
    return %return_values;
}

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

