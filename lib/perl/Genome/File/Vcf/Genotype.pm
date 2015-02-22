package Genome::File::Vcf::Genotype;

use strict;
use warnings;
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

    my @allele_indexes = split(/[\/\|]/, $self->{_gt});
    for my $index (@allele_indexes) {
        unless ($index =~ /[(\d+).]/) {
            confess "Non-numeric allele-index detected: $index";
        }
    }
    $self->{_allele_indexes} = \@allele_indexes;

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
    my $prev;
    for my $index (@{$self->{_allele_indexes}}) {
        unless (defined $prev) {
            $prev = $index;
            next;
        }
        if ($prev != $index) {
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
    return grep {$_ =~ /\./} @{$self->{_allele_indexes}};
}

sub ploidy {
    my $self = shift;
    return scalar @{$self->{_allele_indexes}};
}

sub has_wildtype {
    my $self = shift;
    if ($self->is_missing) {
        return 0;
    }
    for my $index (@{$self->{_allele_indexes}}) {
        if ($index == 0) {
            return 1;
        }
    }
    return 0;
}

sub has_variant {
    my $self = shift;
    if ($self->is_missing) {
        return 0;
    }
    for my $index (@{$self->{_allele_indexes}}) {
        if ($index != 0) {
            return 1;
        }
    }
    return 0;
}

sub get_allele_indexes {
    my $self = shift;
    return @{$self->{_allele_indexes}};
}

sub get_alleles {
    my $self = shift;

    return () if $self->is_missing;

    my @alleles;
    for my $index (@{$self->{_allele_indexes}}) {
        if ($index == 0) {
            push @alleles, $self->{_ref_allele};
        } else {
            push @alleles, $self->{_alt_alleles}->[$index - 1];
        }
    }
    return @alleles;
}

1;

