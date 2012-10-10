package Genome::File::Vep::Entry;

use Carp qw/confess/;
use Genome;
use Storable qw/dclone/;

use strict;
use warnings;

# column offsets in a vep file
our @FIELD_NAMES = qw(
    uploaded_variation
    location
    allele
    gene
    feature
    feature_type
    consequence
    cdna_position
    cds_position
    protein_position
    amino_acids
    codons
    existing_variation
    extra
);

sub new {
    my ($class, $line) = @_;
    chomp $line;

    return bless _parse($line), $class;
}

sub _parse {
    my $line = shift;
    my $self = {};
    @{$self}{@FIELD_NAMES} = split("\t", $line);

    # create new fields for chromosome and position, and parse extra into hash
    ($self->{chrom}, $self->{position}) = split(":", $self->{location});
    $self->{_extra_order} = [
        grep {$_ ne '-' }
        map { s/=.*//; $_ } split(";", $self->{extra}) ];
    
    $self->{extra} = {
        map { split("=") } grep {$_ ne '-'} split(";", $self->{extra})
    };

    return $self;
}

sub to_string {
    my $self = shift;
    my @fields = @{$self}{@FIELD_NAMES};
    $fields[-1] = join(";", map {"$_=$self->{extra}{$_}"} @{$self->{_extra_order}}) || '-';
    return join("\t", @fields);
}

sub set_extra_field {
    my ($self, $key, $value) = @_;
    if (defined $value) {
        push(@{$self->{_extra_order}}, $key) unless grep {$_ eq $key} @{$self->{_extra_order}};
        $self->{extra}{$key} = $value;
    } else {
        @{$self->{_extra_order}} = grep {$_ ne $key} @{$self->{_extra_order}};
        delete $self->{extra}{$key};
    }
}

sub clone {
    my $self = shift;
    return dclone($self);
}

1;
