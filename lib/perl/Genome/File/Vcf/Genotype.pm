package Genome::File::Vcf::Genotype;

use strict;
use warnings;
use Genome;

sub new {
    my ($class, $ref_allele, $alt_alleles, $gt) = @_;
    my $self = {
        _ref_allele => $ref_allele,
        _alt_alleles => $alt_alleles,
        _gt => $gt,
    };
    bless $self, $class;
    $self->_parse;
    return $self;
}

sub _parse {
    my $self = shift;
    my @alleles = split("/", $self->{_gt});
    my $num_alleles = scalar @alleles;
    for my $index (0..$num_alleles-1) {
        unless ($alleles[$index] =~ /[0-9]+/) {
            delete $alleles[$index];
        }
    }
    $self->{_alleles} = \@alleles;
}

sub is_homozygous {
    my $self = shift;
    if ($self->is_missing) {
        return 0;
    }
    my $is_homozygous = 0;
    my $prevAllele;
    for my $allele (@{$self->{_alleles}}) {
        unless (defined $prevAllele) {
            $prevAllele = $allele;
            next;
        }
        if ($prevAllele != $allele) {
            return 0;
        }
    }
    return 1;
}

sub is_heterozygous {
    my $self = shift;    
    if ($self->is_missing) {
        return 0;
    }
    return !$self->is_homozygous;
}

sub is_missing {
    my $self = shift;
    my $num_alleles = scalar @{$self->{_alleles}};
    if ($num_alleles <= 0) {
        return 1;
    }
    return 0;
}

sub is_wildtype {
    my $self = shift;
    if ($self->is_missing) {
        return 0;
    }for my $allele (@{$self->{_alleles}}) {
        if ($allele == 0) {
            return 1;
        }
    }
    return 0;
}

sub is_variant {
    my $self = shift;
    if ($self->is_missing) {
        return 0;
    }for my $allele (@{$self->{_alleles}}) {
        if ($allele != 0) {
            return 1;
        }
    }
    return 0;
}

sub get_allele_by_index {
    my $self = shift;
    my $index = shift;
    return $self->{_alleles}->[$index];
}

1;

