package Genome::File::Vcf::Differ;

use strict;
use warnings;

use Genome;
use Genome::File::Vcf::Reader;
use Genome::File::Vcf::EntryDiff;
use Genome::File::Vcf::HeaderDiff;
use Memoize;
use Carp qw/confess/;
use Set::Scalar;
use Genome::File::Vcf::HeaderDiff;

=head1 NAME

Genome::File::Vcf::Differ - A class for diffing vcf files.

=head1 SYNOPSIS

use Genome::File::Vcf::Differ;

my $differ = new Genome::File::Vcf::Differ("A.vcf", "B.vcf.gz")
while (my $diff = $differ->diff) {
    # $diff is either a header diff or an entry diff
    $diff->print;
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
            @columns = $b->{header}->all_columns;
        } elsif (!defined($b)) { # a has MORE lines than b
            @columns = $a->{header}->all_columns;
        } else {
            @columns = $self->entries_diff($a, $b);
        }

        if (@columns) {
            return Genome::File::Vcf::EntryDiff->new(
                $self->{_a} => $a,
                $self->{_b} => $b,
                \@columns);
        }
    }
    return;
}

sub entries_diff {
    my ($self, $a, $b) = @_;

    my %a_hash = %{$a->to_hashref};
    my %b_hash = %{$b->to_hashref};

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

sub header_a {
    my $self = shift;
    return $self->{_a_reader}->header;
}

sub header_b {
    my $self = shift;
    return $self->{_b_reader}->header;
}

sub header {
    my $self = shift;

    my ($lines_a, $samples_a) = _lines_and_samples_from_header($self->header_a);
    my ($lines_b, $samples_b) = _lines_and_samples_from_header($self->header_b);

    my $line_diffs_a = $lines_a->difference($lines_b);
    my $line_diffs_b = $lines_b->difference($lines_a);

    my $sample_diffs_a = $samples_a->difference($samples_b);
    my $sample_diffs_b = $samples_b->difference($samples_a);

    if (all_empty($line_diffs_a, $line_diffs_b, $sample_diffs_a, $sample_diffs_b)) {
        return;
    } else {
        return Genome::File::Vcf::HeaderDiff->new(
            $self->{_a}, $line_diffs_a, $sample_diffs_a,
            $self->{_b}, $line_diffs_b, $sample_diffs_b,
        );
    }
}

sub _lines_and_samples_from_header {
    my $header = shift;

    my $lines = Set::Scalar->new();
    $lines->insert($header->_info_lines);
    $lines->insert($header->_format_lines);
    $lines->insert($header->_filter_lines);

    my $samples = Set::Scalar->new($header->sample_names);
    return ($lines, $samples);
}

sub all_empty {
    for my $set (@_) {
        if ($set->size() > 0) {
            return 0;
        }
    }
    return 1;
}

1;
