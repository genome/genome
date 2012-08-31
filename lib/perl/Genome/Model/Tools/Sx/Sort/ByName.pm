package Genome::Model::Tools::Sx::Sort::ByName;

use strict;
use warnings;
no warnings 'uninitialized';

use Genome;

class Genome::Model::Tools::Sx::Sort::ByName {
    is  => 'Genome::Model::Tools::Sx::Sort::Base',
};

sub help_brief {
    return 'Sort sequences by name';
}

sub _sort_params {
    return (qw/ -k1 /);
}

sub _flatten {
    # tab separated seqs on one line
    return join("\t", map { my $seq = $_; map { $seq->{$_} } (qw/ id desc seq qual /) } @_);
}

sub _inflate {
    chomp $_[0];
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

