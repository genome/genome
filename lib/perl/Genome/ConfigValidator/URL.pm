package Genome::ConfigValidator::URL;

use strict;
use warnings;

require URI;

use Mouse;

with qw(Genome::ConfigValidatorBase);

sub check {
    my ($self, $value) = @_;
    my $uri = URI->new($value);
    return $uri->as_string eq $value
        && $uri->scheme;
}

sub message {
    return 'an URL';
}

1;
