package Genome::ConfigValidator::URI;

use strict;
use warnings;

require URI;

use Mouse;

with qw(Genome::ConfigValidatorBase);

sub check {
    my ($self, $value) = @_;
    my $uri = URI->new($value);
    return $uri->as_string eq $value;
}

sub message {
    return 'an URI';
}

1;
