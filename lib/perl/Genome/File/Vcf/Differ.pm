package Genome::File::Vcf::Differ;

use strict;
use warnings;

use Genome;
use Genome::File::Vcf::Reader;
use Carp qw/confess/;
use Set::Scalar;

=head1 NAME

Genome::File::Vcf::Differ - A class for diffing vcf files.

=head1 SYNOPSIS

my $differ = new Genome::File::Vcf::Differ("A.vcf", "B.vcf.gz")

#   keys = file_names
#   values = list of lines that differ
my %header_differences = $differ->header;


while (my ($A_entry, $B_entry, @columns) = $differ->next) {
    # ...
}

=cut

sub new {
    my ($class, $a, $b) = @_;

    my $a_reader = Genome::File::Vcf::Reader->new($a);
    my $b_reader = Genome::File::Vcf::Reader->new($b);

    my @base_columns = qw(CHROM POS ID REF ALT QUAL FILTER INFO FORMAT);
    my @a_columns = (@base_columns, $a_reader->header->sample_names);
    my @b_columns = (@base_columns, $b_reader->header->sample_names);

    my $self = {
        _a => $a,
        _b => $b,
        _a_reader => $a_reader,
        _b_reader => $b_reader,
        _a_columns => \@a_columns,
        _b_columns => \@b_columns,
    };

    bless $self, $class;
    return $self;
}

sub diff {
    my $self = shift;

    if ($self->header) {
        return $self->header;
    }
    return $self->next;
}

sub next {
    my $self = shift;

    while (my ($a, $b) = $self->paired_next) {
        my @columns;

        if (!defined($a)) { # a has FEWER lines than b
            @columns = @{$self->{_b_columns}};
        } elsif (!defined($b)) { # a has MORE lines than b
            @columns = @{$self->{_a_columns}};
        } else {
            @columns = $self->entries_diff($a, $b);
        }

        if (@columns) {
            return ($a, $b, @columns);
        }
    }
    return;
}

sub entries_diff {
    my ($self, $a, $b) = @_;

    my @a_keys = @{$self->{_a_columns}};
    my @a_values = split(/\t/, $a->to_string);
    my %a_hash;
    @a_hash{@a_keys} = @a_values;

    my @b_keys = @{$self->{_b_columns}};
    my @b_values = split(/\t/, $b->to_string);
    my %b_hash;
    @b_hash{@b_keys} = @b_values;

    my @columns;
    for my $key (@a_keys) {
        if ($a_hash{$key} ne $b_hash{$key}) {
            push @columns, $key;
        }
    }
    return @columns;
}

sub paired_next {
    my $self = shift;

    my $next_a = $self->{_a_reader}->next;
    my $next_b = $self->{_b_reader}->next;

    if (!defined($next_a) and !defined($next_b)) {
        return;
    } else {
        return ($next_a, $next_b);
    }
}

sub header {
    my $self = shift;

    my $ah = $self->{_a_reader}->header;
    my $bh = $self->{_b_reader}->header;

    my $a_lineset = _lineset_from_header($ah);
    my $b_lineset = _lineset_from_header($bh);

    my $a_diff = $a_lineset->difference($b_lineset);
    my $b_diff = $b_lineset->difference($a_lineset);

    if (!$a_diff->size && !$b_diff->size) {
        return;
    } else {
        my @a_diffs = sort($a_diff->members());
        my @b_diffs = sort($b_diff->members());

        return (
            $self->{_a} => \@a_diffs,
            $self->{_b} => \@b_diffs,
        );
    }
}

sub _lineset_from_header {
    my $header = shift;
    my $lineset = Set::Scalar->new();
    $lineset->insert($header->_info_lines);
    $lineset->insert($header->_format_lines);
    $lineset->insert($header->_filter_lines);
    $lineset->insert($header->sample_names);
    return $lineset;
}

1;
