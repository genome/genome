#!/usr/bin/env genome-perl

use strict;
use warnings;

use above "Genome";

use Test::More tests => 2;

use Genome::Utility::Email;

{
    my $addr = Genome::Utility::Email::construct_address();
    my $expected = sprintf('%s@%s', $ENV{USER}, $ENV{GENOME_EMAIL_DOMAIN});
    is($addr, $expected, 'construct_address with no params');
}

{
    my $name = 'bob';
    $name++ while ($name eq $ENV{USER});
    my $addr = Genome::Utility::Email::construct_address($name);
    my $expected = sprintf('%s@%s', $name, $ENV{GENOME_EMAIL_DOMAIN});
    is($addr, $expected, 'construct_address given a name');
}
