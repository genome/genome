package Genome::File::BedPe::Header;

use strict;
use warnings;

use Exporter;
use base qw(Exporter);

our @REQUIRED_FIELDS = qw{
    chrom1
    start1
    end1
    chrom2
    start2
    end2
};

our @OPTIONAL_FIELDS = qw{
    name
    score
    strand1
    strand2
};

our @ALL_FIELDS = (@REQUIRED_FIELDS, @OPTIONAL_FIELDS);

our %FIELD_INDICES = map {$ALL_FIELDS[$_] => $_} 0..$#ALL_FIELDS;

our @EXPORT_OK = qw(@REQUIRED_FIELDS @EXTRA_FIELDS @ALL_FIELDS %FIELD_INDICES);

sub new {
    my ($class, $lines) = @_;

    $lines = [] unless defined $lines;

    my $self = {
        lines => $lines,
        custom_fields => [],
        _custom_field_idx => {},
    };

    bless $self, $class;
    return $self;
}

sub make_header_line {
    my $self = shift;
    die "make_header_line called on header that already has data!" if @{$self->{lines}};
    push @{$self->{lines}},
        join("\t", @ALL_FIELDS, @{$self->{custom_fields}});
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

sub to_string {
    my $self = shift;

    return join("\n", @{$self->{lines}});
}

1;
