package Genome::File::BedPe::Header;

use strict;
use warnings;

sub new {
    my ($class, $lines) = @_;

    my $self = {
        lines => $lines,
        custom_fields => [],
        _custom_field_idx => {},
    };

    bless $self, $class;
    return $self;
}

sub _build_custom_field_idx {
    my $self = shift;
    $self->{_custom_field_idx} = {
        map { $self->{custom_fields}[$_] => $_ } 0..$#{$self->{custom_fields}}
    };
}

sub set_custom_fields {
    my ($self, @fields) = @_;
    $self->{custom_fields} = \@fields;
    $self->_build_custom_field_idx;
}

sub guess_custom_fields {
    my $self = shift;
    if (@{$self->{lines}}) {
        my $last_line = $self->{lines}[-1];
        my @fields = split("\t", $last_line);
        if ($#fields >= 10) {
            $self->{custom_fields} = [@fields[10..$#fields]];
            $self->_build_custom_field_idx;
            return 1;
        }
    }
    return 1;
}

sub custom_field_index {
    my ($self, $field_name) = @_;
    if (exists $self->{_custom_field_idx}{$field_name}) {
        return $self->{_custom_field_idx}{$field_name};
    }
    return;
}

1;
