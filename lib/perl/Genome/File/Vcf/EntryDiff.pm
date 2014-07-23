package Genome::File::Vcf::EntryDiff;

use strict;
use warnings FATAL => 'all';
use Params::Validate qw(validate validate_pos :types);

=head1 NAME

Genome::File::Vcf::EntryDiff - Representation of differences between
two VCF entries.

=head1 SYNOPSIS

    use Genome::File::Vcf::Differ;

    my $differ = Genome::File::Vcf::Differ->new('filea.vcf', 'fileb.vcf');
    my $diff = $differ->next;
    if ($diff) {
       $diff->print;
    }

=head1 DESCRIPTION

    The EntryDiff print method can display one of three categories:

    1) The differences between columns when both entries represent the same position
    2) The different chromosomes and positions between entries
    3) One file contains more lines than the other

=cut

sub new {
    my ($class, $a_path, $entry_a, $b_path, $entry_b, $columns) = validate_pos(@_,
        {type => SCALAR},
        {type => SCALAR},
        {type => OBJECT | UNDEF},
        {type => SCALAR},
        {type => OBJECT | UNDEF},
        {type => ARRAYREF},
    );
    my $self = {
        _entry_a => $entry_a,
        _a => $a_path,
        _entry_b => $entry_b,
        _b => $b_path,
        _columns => $columns
    };
 
    bless $self, $class;
    return $self;
}

sub print {
    my $self = shift;

    print $self->to_string . "\n";
}

sub to_string {
    my $self = shift;
    my $string_rep;

    if ($self->_both_entries_defined) {
        $string_rep = $self->_explain_entry_diffs;
    } else {
        $string_rep = $self->_explain_entry_not_defined;
    }
    
    return $string_rep;
}

sub _both_entries_defined {
    my $self = shift;

    return ( defined($self->{_entry_a}) && defined($self->{_entry_b}) );
}

sub _explain_entry_diffs {
    my $self = shift;
    my $entry_diffs;

    if ($self->_both_same_chrom_pos) {
       $entry_diffs = $self->_explain_column_diffs;
    } else {
       $entry_diffs =$self->_explain_position_diffs;
    }

    return $entry_diffs;
}

sub _both_same_chrom_pos {
    my $self = shift;

    return  ( ($self->{_entry_a}->{chrom} eq $self->{_entry_b}->{chrom}) &&
              ($self->{_entry_a}->{position} ==  $self->{_entry_b}->{position}) );
}

sub _explain_column_diffs {
    my $self = shift;
    my $column_diffs;

    $column_diffs = sprintf "Entries at position (%s, %d) have differences in columns: %s\n",
        $self->{_entry_a}->{chrom}, $self->{_entry_a}->{position}, join(', ', @{$self->{_columns}});
    $column_diffs .= sprintf "a: %s, b: %s\n", $self->{_a}, $self->{_b};

    for my $column (@{$self->{_columns}}) {
        $column_diffs .= sprintf "%s:\n    a => %s\n    b => %s\n", $column, 
        $self->{_entry_a}->to_hashref->{$column},
        $self->{_entry_b}->to_hashref->{$column};
    }
    return $column_diffs;
}

sub _explain_position_diffs {
    my $self = shift;
    my $pos_diffs;

    $pos_diffs = sprintf "Entries refer to different chromosomes and positions:\n" .
           "%10s %12d %s\n" .
           "%10s %12d %s\n",
           $self->{_entry_a}->{chrom}, $self->{_entry_a}->{position}, $self->{_a},
           $self->{_entry_b}->{chrom}, $self->{_entry_b}->{position}, $self->{_b};

    return $pos_diffs;
}

sub _explain_entry_not_defined {
    my $self = shift;
    my $short_entry;

    if ( !defined($self->{_entry_a}) ) {
        $short_entry = sprintf "%s has fewer lines than %s\n", $self->{_a}, $self->{_b};
    } else {
        $short_entry = sprintf "%s has fewer lines than %s\n", $self->{_b}, $self->{_a};
    }
    
    return $short_entry;
}

1;
