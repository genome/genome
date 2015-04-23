package Genome::ConfigValidator::URI;

use strict;
use warnings;

require URI;

sub validate {
    my $value = shift;
    my $uri = URI->new($value);
    return if $uri->as_string eq $value;
    return 'a URI';
}

1;

