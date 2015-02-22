#!/usr/bin/env genome-perl

use strict;
use warnings FATAL => 'all';

use Test::More;
use above 'Genome';

my $pkg = 'Genome::File::Location';
use_ok($pkg) || die;

my $a = $pkg->new(1, 1);
my $b = $pkg->new('X', 1);

expect_less($a, $b);

$a->set_equal_to($b);
expect_equal($a, $b);

$b->{pos} = 2;
expect_less($a, $b);

done_testing();

sub expect_equal {
    my ($a, $b) = @_;

    is($a->is_greater_than($b), '',
        sprintf("%s >  %s is false", $a->to_string, $b->to_string));
    is($a->is_greater_than_or_equal_to($b), 1,
        sprintf("%s >= %s is true", $a->to_string, $b->to_string));
    is($a->is_less_than($b), '',
        sprintf("%s <  %s is false", $a->to_string, $b->to_string));
    is($a->is_less_than_or_equal_to($b), 1,
        sprintf("%s <= %s is true", $a->to_string, $b->to_string));
    is($a->is_equal_to($b), 1,
        sprintf("%s == %s is true", $a->to_string, $b->to_string));
}

sub expect_less {
    my ($a, $b) = @_;

    is($a->is_greater_than($b), '',
        sprintf("%s >  %s is false", $a->to_string, $b->to_string));
    is($a->is_greater_than_or_equal_to($b), '',
        sprintf("%s >= %s is false", $a->to_string, $b->to_string));
    is($a->is_less_than($b), 1,
        sprintf("%s <  %s is true", $a->to_string, $b->to_string));
    is($a->is_less_than_or_equal_to($b), 1,
        sprintf("%s <= %s is true", $a->to_string, $b->to_string));
    is($a->is_equal_to($b), '',
        sprintf("%s == %s is false", $a->to_string, $b->to_string));
}

