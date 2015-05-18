package Genome::ConfigValidator::LSFQueue;

use strict;
use warnings;

use Mouse;

with qw(Genome::ConfigValidatorBase);

sub check {
    my ($self, $value) = @_;
    system(qq(bqueues $value 1> /dev/null 2>&1));
    return $? == 0;
}

sub message {
    return 'an LSF queue';
}

1;
