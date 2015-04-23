package Genome::ConfigValidator::LSFQueue;

use strict;
use warnings;

require URI;

sub validate {
    my $value = shift;
    system(qq(bqueues $value 1> /dev/null 2>&1));
    return if $? == 0;
    return 'a valid LSF queue';
}

1;
