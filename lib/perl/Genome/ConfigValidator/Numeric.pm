package Genome::ConfigValidator::Numeric;

use strict;
use warnings;

use Scalar::Util qw(looks_like_number);

sub validate {
    my $value = shift;
    return if looks_like_number($value);
    return 'numeric';
}

1;
