package Genome::ConfigValidator::Positive;

use strict;
use warnings;

use Scalar::Util qw(looks_like_number);

sub validate {
    my $value = shift;
    return if ( looks_like_number($value) && $value > 0 );
    return 'positive';
}

1;
