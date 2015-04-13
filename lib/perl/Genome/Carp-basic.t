#!/usr/bin/env genome-perl

use strict;
use warnings;

use Test::More tests => 9;
use Test::Fatal qw(exception);

use_ok('Genome::Carp');

my @die_methods = qw(croakf confessf dief);
for my $die_method (@die_methods) {
    my $method = \&{ "Genome::Carp::$die_method" };
    my $exception = exception { $method->('hello %s', 'joe') };
    like($exception, qr/^hello joe/, "$die_method did die");
}

my @warn_methods = qw(carpf cluckf warnf);
for my $warn_method (@warn_methods) {
    my $method = \&{ "Genome::Carp::$warn_method" };
    my $warning;
    local $SIG{__WARN__} = sub { $warning = $_[0] };
    $method->('hello %s', 'joe');
    like($warning, qr/^hello joe/, "$warn_method did warn");
}

my @string_methods = qw(longmessf shortmessf);
for my $string_method (@string_methods) {
    my $method = \&{ "Genome::Carp::$string_method" };
    my $got = $method->('hello %s', 'joe');
    like($got, qr/^hello joe/, "$string_method did match");
}
