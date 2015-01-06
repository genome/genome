package Genome::File::Location;

use strict;
use warnings FATAL => 'all';
use Data::Dump qw(pp);
use Sort::strverscmp;

sub new {
    my ($class, $chrom, $pos) = @_;

    my $self = {
        chrom => $chrom,
        pos => $pos,
    };

    bless $self, $class;
    return $self;
}

sub set_equal_to {
    my ($self, $other) = @_;
    $self->{chrom} = $other->{chrom};
    $self->{pos} = $other->{pos};
    return;
}

sub cmp_with {
    my ($self, $other) = @_;

    my $chrom_cmp = strverscmp($self->{chrom}, $other->{chrom});
    if ($chrom_cmp == 0) {
        return strverscmp($self->{pos}, $other->{pos});
    } else {
        return $chrom_cmp;
    }
}

sub is_greater_than {
    my ($self, $other) = @_;
    return $self->cmp_with($other) == 1;
}

sub is_greater_than_or_equal_to {
    my ($self, $other) = @_;
    return $self->cmp_with($other) > -1;
}

sub is_less_than {
    my ($self, $other) = @_;
    return $self->cmp_with($other) == -1;
}

sub is_less_than_or_equal_to {
    my ($self, $other) = @_;
    return $self->cmp_with($other) < 1;
}

sub is_equal_to {
    my ($self, $other) = @_;
    return $self->cmp_with($other) == 0;
}

sub to_string {
    my $self = shift;
    return sprintf("(%s, %s)", pp($self->{chrom}), pp($self->{pos}));
}

1;
