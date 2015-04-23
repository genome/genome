package Genome::ConfigValidator::URL;

use strict;
use warnings;

require URI;

sub validate {
    my $value = shift;
    my $uri = URI->new($value);
    return if $uri->scheme;
    return 'a URL';
}

1;
