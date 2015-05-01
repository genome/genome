package Genome::ConfigValidator::Defined;

use strict;
use warnings;

use Mouse;

with qw(Genome::ConfigValidatorBase);

sub check {
    my ($self, $value) = @_;
    return defined($value)
}

sub message {
    return 'defined';
}

1;
