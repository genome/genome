package Genome::File::Vcf::Differ;

use strict;
use warnings;

use Genome;
use Genome::File::Vcf::Reader;
use Memoize;
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

    my $self = {
        _a => $a,
        _b => $b,
        _a_reader => $a_reader,
        _b_reader => $b_reader,
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
            @columns = $self->columns_for_header($b->{header});
        } elsif (!defined($b)) { # a has MORE lines than b
            @columns = $self->columns_for_header($a->{header});
        } else {
            @columns = $self->entries_diff($a, $b);
        }

        if (@columns) {
            return ($a, $b, @columns);
        }
    }
    return;
}

sub columns_for_header {
    my ($self, $header) = @_;

    my @base_columns = qw(CHROM POS ID REF ALT QUAL FILTER INFO FORMAT);
    my @columns = (@base_columns, $header->sample_names);
    return @columns;
}
Memoize::memoize('columns_for_header');

sub hash_repr {
    my ($self, $entry) = @_;

    my @keys = $self->columns_for_header($entry->{header});
    my @values = split(/\t/, $entry->to_string);

    my %hash;
    @hash{@keys} = @values;
    return %hash;
}

sub entries_diff {
    my ($self, $a, $b) = @_;

    my %a_hash = $self->hash_repr($a);
    my %b_hash = $self->hash_repr($b);

    my @columns;
    for my $key (keys %a_hash) {
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
