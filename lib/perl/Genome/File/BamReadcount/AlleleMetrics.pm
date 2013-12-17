package Genome::File::BamReadcount::AlleleMetrics;

use Genome;

use strict;
use warnings;

my @properties = qw/
_allele
_count
_avg_mapq
_avg_bq
_avg_se_mapq
_num_plus_strand
_num_minus_strand
_avg_pos_as_fraction
_avg_num_mismatches_as_fraction
_avg_sum_mismatch_qualities
_num_q2_containing_reads
_avg_distance_to_q2_start_in_q2_reads
_avg_clipped_length
_avg_distance_to_effective_3p_end
/;

sub new {
    my ($class, $string) = @_;
    my %property_hash;
    $property_hash{_source_string} = $string;
    @property_hash{@properties} = split /:/, $string;
    for my $property (@properties) {
        unless(defined $property_hash{$property}) {
            die "Invalid entry\n";
        }
    }
    my $self = \%property_hash;
    bless $self, $class;
    return $self;
}

sub allele {
    my ($self) = @_;
    return $self->{_allele};
}

sub is_indel {
    my ($self) = @_;
    return ($self->is_del or $self->is_ins);
}

sub is_del {
    my ($self) = @_;
    return $self->{_allele} =~ /^[-]/;
}

sub is_ins {
    my ($self) = @_;
    return $self->{_allele} =~ /^[+]/;
}

sub count {
    my ($self) = @_;
    return $self->{_count};
}

sub avg_mapq {
    my ($self) = @_;
    return $self->{_avg_mapq};
}

sub avg_bq {
    my ($self) = @_;
    return $self->{_avg_bq};
}

sub avg_se_mapq {
    my ($self) = @_;
    return $self->{_avg_se_mapq};
}

sub num_plus_strand {
    my ($self) = @_;
    return $self->{_num_plus_strand};
}

sub num_minus_strand {
    my ($self) = @_;
    return $self->{_num_minus_strand};
}

sub avg_pos_as_fraction {
    my ($self) = @_;
    return $self->{_avg_pos_as_fraction};
}

sub avg_num_mismatches_as_fraction {
    my ($self) = @_;
    return $self->{_avg_num_mismatches_as_fraction};
}

sub avg_sum_mismatch_qualities {
    my ($self) = @_;
    return $self->{_avg_sum_mismatch_qualities};
}

sub num_q2_containing_reads {
    my ($self) = @_;
    return $self->{_num_q2_containing_reads};
}

sub avg_distance_to_q2_start_in_q2_reads {
    my ($self) = @_;
    return $self->{_avg_distance_to_q2_start_in_q2_reads};
}

sub avg_clipped_length {
    my ($self) = @_;
    return $self->{_avg_clipped_length};
}

sub avg_distance_to_effective_3p_end {
    my ($self) = @_;
    return $self->{_avg_distance_to_effective_3p_end};
}

1;
