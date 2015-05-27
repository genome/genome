package Genome::ConfigValidator::Numeric;

use strict;
use warnings;

require Scalar::Util;

use Mouse;

with qw(Genome::ConfigValidatorBase);

sub check {
    my ($self, $value) = @_;
    return Scalar::Util::looks_like_number($value);
}

sub message {
    return 'a number';
}

1;
