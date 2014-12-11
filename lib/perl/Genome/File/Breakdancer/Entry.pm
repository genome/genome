package Genome::File::Breakdancer::Entry;

use strict;
use warnings;

my @FIELD_NAMES = qw(
    chr1 pos1 orientation1
    chr2 pos2 orientation2
    type size score
    num_reads
    num_reads_lib
    );


sub new {
    my ($class, $header, $line) = @_;
    my $self = {
        header => $header,
        _line => $line,
        _lib_order => [],
        extra => [],
    };
    bless $self, $class;
    $self->_parse;
    return $self;
}

sub _parse_orientation {
    my $ori = shift;
    my ($pos) = $ori =~ /([0-9]+)\+/;
    my ($neg) = $ori =~ /([0-9]+)\-/;
    return {fwd => $pos || 0, rev => $neg || 0};
}

sub _parse_lib_counts_cn {
    my $str = shift;
    my ($lib, $count, $cn) = split(/[,|]/, $str);
    return $lib => {count => $count, copy_number => $cn};
}

sub _format_lib_counts_cn {
    my $self = shift;

    my $libcounts = $self->{num_reads_lib};
    return join(":", map({
            sprintf("%s|%d,%s", $_, $libcounts->{$_}{count},
                $libcounts->{$_}{copy_number});
            } @{$self->{_lib_order}})
        );
}

sub _parse {
    my ($self) = @_;

    my @fields = split("\t", $self->{_line});
    @{$self}{@FIELD_NAMES} = @fields[0..$#FIELD_NAMES];

    my @lib_read_counts = map {_parse_lib_counts_cn($_)} split(":", $self->{num_reads_lib});
    @{$self->{_lib_order}} = @lib_read_counts[grep {$_ % 2 == 0} 0..$#lib_read_counts];

    $self->{num_reads_lib} = {@lib_read_counts};
    $self->{$_} = _parse_orientation($self->{$_}) for qw/orientation1 orientation2/;
}

sub lib_read_count {
    my ($self, $lib) = @_;
    if (exists $self->{num_reads_lib}{$lib}) {
        return $self->{num_reads_lib}{$lib}{count};
    }
    else {
        return 0;
    }
}

sub to_string {
    my $self = shift;

    return join("\t",
        $self->{chr1}, $self->{pos1},
        sprintf("%d+%d-", $self->{orientation1}{fwd}, $self->{orientation1}{rev}),
        $self->{chr2}, $self->{pos2},
        sprintf("%d+%d-", $self->{orientation2}{fwd}, $self->{orientation2}{rev}),
        $self->{type}, $self->{size}, $self->{score}, $self->{num_reads},
        $self->_format_lib_counts_cn,
        @{$self->{extra}});
}

1;
