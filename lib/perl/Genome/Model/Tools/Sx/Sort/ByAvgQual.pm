package Genome::Model::Tools::Sx::Sort::ByAvgQual;

use strict;
use warnings;
no warnings 'uninitialized'; # desc may be undefined for sequences

use Genome;

require List::Util;

class Genome::Model::Tools::Sx::Sort::ByAvgQual {
    is  => 'Genome::Model::Tools::Sx::Sort::Base',
};

sub help_brief {
    return 'Sort sequences by average quality';
}

sub _sort_params {
    return (qw/ -nr -k1 /);
}

sub _flatten {
    # tab separated, avg qual first, then seqs
    return join(
        "\t",
        List::Util::sum( map { Genome::Model::Tools::Sx::Functions->calculate_average_quality($_->{qual}) } @_ ),
        map { my $seq = $_; map { $seq->{$_} } (qw/ id desc seq qual /) } @_,
    );
}

sub _inflate {
    chomp $_[0];
    $_[0] =~ s/\d+\t//; # rm avg qual
    my @tokens = split("\t", $_[0]);
    my @seqs;
    for ( my $i = 0; $i < $#tokens; $i += 4 ) {
        my %seq = (
            id => $tokens[$i],
            desc => $tokens[$i + 1],
            seq => $tokens[$i + 2],
            qual => $tokens[$i + 3],
        );
        push @seqs, \%seq;
    }
    return \@seqs
}

1;

