use strict;
use warnings;
use Test::More;

use above "Genome";

use Genome::Model::Tools::Analysis::Helpers qw(
    byBamOrder
);

my @unsorted = (
    "X\t44",
    "Y\t50",
    "1\t100",
    "1\t1",
);

my @expected = (
    "1\t1",
    "1\t100",
    "X\t44",
    "Y\t50",
);

my @sorted = sort byBamOrder @unsorted;
is_deeply(\@sorted, \@expected, 'Sorted byBamOrder works as expected');

done_testing();

