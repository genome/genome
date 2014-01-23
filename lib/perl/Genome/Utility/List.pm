package Genome::Utility::List;

use strict;
use warnings;

use List::Util 'first';

our @ISA = qw(Exporter);
our @EXPORT_OK = qw(
        in join_with_single_slash
);

sub in {
    my $element = shift;
    my $find = first {$element eq $_} @_;
    return defined($find);
}

sub join_with_single_slash {
    List::Util::reduce
        {   $a =~ s/\/$//;
            $b =~ s/^\///;
            return join('/', $a, $b);
        }
        @_;
}

1;
