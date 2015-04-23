package Genome::ConfigValidator::Port;

use strict;
use warnings;

use Scalar::Util qw(looks_like_number);

sub validate {
    my $value = shift;
    if (looks_like_number($value)
        && $value <= 65535
        && $value >= 1
    ) {
        return;
    }
    return 'a Port';
}

1;

