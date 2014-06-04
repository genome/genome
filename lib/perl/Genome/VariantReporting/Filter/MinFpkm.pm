package Genome::VariantReporting::Filter::MinFpkm;

use strict;
use warnings;
use Genome;
use Scalar::Util qw(looks_like_number);

class Genome::VariantReporting::Filter::MinFpkm {
    is => ['Genome::VariantReporting::Filter::Base', 'Genome::VariantReporting::WithSampleName'],
    has => [
        min_fpkm => {
            is => "Number",
            doc => "The minimum FPKM value to pass",
        },
    ],
};

sub name {
    return 'min-fpkm';
}

sub __errors__ {
    my $self = shift;
    my @errors = $self->SUPER::__errors__;

    my $min_fpkm = $self->min_fpkm;
    unless ( ( looks_like_number($min_fpkm) ) and ($min_fpkm >= 0) ){
        push @errors, UR::Object::Tag->create(
            type => 'error',
            properties => [qw/ min_fpkm /],
            desc => "Value given for min_fpkm ($min_fpkm) is not a positive number!",
        );
    }

    return @errors;
}

sub filter_entry {
    my $self = shift;
    my $entry = shift;

    my %return_values;
    my %fpkm_for_genotype_allele = $self->fpkm_for_genotype_allele($entry);
    for my $alt_allele ( @{$entry->{alternate_alleles}} ) {
        $return_values{$alt_allele} = $self->value_passes($fpkm_for_genotype_allele{$alt_allele});
    }

    return %return_values;
}

sub value_passes {
    my ($self, $value) = @_;
    if ( (not defined $value ) or ($value eq '.') or ($value < $self->min_fpkm) ) {
        return 0;
    } else {
        return 1;
    }
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

