package Genome::Utility::List;

use strict;
use warnings;

use List::Util 'first';

our @ISA = qw(Exporter);
our @EXPORT_OK = qw(
        in
);

sub in {
    my $element = shift;
    my $find = first {$element eq $_} @_;
    return defined($find);
}
