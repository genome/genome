package Genome::ConfigValidator::Port;

use strict;
use warnings;

require Scalar::Util;

use Mouse;

with qw(Genome::ConfigValidatorBase);

sub check {
    my ($self, $value) = @_;
    return Scalar::Util::looks_like_number($value)
        && $value <= 65535
        && $value >= 1;
}

sub message {
    return 'a port';
}

1;
