package Genome::File::Vcf::Genotype;

use strict;
use warnings;
use Genome;
use Carp qw/confess/;

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
    confess "Attempted to parse null genotype" unless $self->{_gt};

    my @alleles = grep { /[0-9]+/ } split(qr{/|}, $self->{_gt});
    
    $self->{_alleles} = \@alleles;

    $self->{_is_phased} = _is_phased($self->{_gt});
}

sub _is_phased {
    my $gt = shift;
    if ($gt =~ /\|/) {
        return 1;
    }
    return 0;
}

sub is_phased {
    my $self = shift;
    return $self->{_is_phased};
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

sub has_wildtype {
    my $self = shift;
    if ($self->is_missing) {
        return 0;
    }
    for my $allele (@{$self->{_alleles}}) {
        if ($allele == 0) {
            return 1;
        }
    }
    return 0;
}

sub has_variant {
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

sub get_alleles {
    my $self = shift;
    return @{$self->{_alleles}};
}

1;

