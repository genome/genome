package Genome::File::BedPe::Entry;

use Genome;
use Carp qw/confess/;
use List::AllUtils qw/min first/;

use strict;
use warnings;

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

sub new {
    my ($class, $header, $line) = @_;
    my $self = {
        header => $header,
        _line => $line,
        custom => [],
    };

    bless $self, $class;
    $self->_parse;
    return $self;
}

sub _parse {
    my $self = shift;
    my @fields = split("\t", $self->{_line});
    if (scalar @fields < scalar @REQUIRED_FIELDS) {
        confess "Too few fields in bedpe record (" . scalar @fields . ")";
    }

    my $last_field = min($#fields, $#ALL_FIELDS);

    @{$self}{@ALL_FIELDS[0..$last_field]} = @fields[0..$last_field];
    if ($#fields > $last_field) {
        $self->{custom} = [ @fields[$last_field + 1 .. $#fields] ];
    }
}

sub validate {
    my $self = shift;
    # numeric fields, (-1 means not known)
    my @numeric_fields = qw/start1 start2 end1 end2/;

    for my $f (@numeric_fields) {
        next if $f eq '-1';
        if ($self->{$f} !~ /^\d+$/) {
            confess "in entry $self->{_line}: $f is not numeric";
        }
    }
}

sub to_string {
    my $self = shift;
    my $last_field = first { !defined $self->{$ALL_FIELDS[$_]} } 0..$#ALL_FIELDS;
    if (!defined $last_field) {
        $last_field = $#ALL_FIELDS;
    }
    else {
        --$last_field;
    }

    return join("\t", @{$self}{@ALL_FIELDS[0..$last_field]}, @{$self->{custom}});
}

1;
