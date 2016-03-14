#!/usr/bin/env genome-perl

use strict;
use warnings;

BEGIN {
    $ENV{UR_DBI_NO_COMMIT} = 1;
}

use above "Genome";

use Test::More tests => 2;

use Genome::Utility::Email;

{
    my $me = Genome::Sys->current_user;
    my $addr = Genome::Utility::Email::construct_address();
    my $expected = $me->email;
    is($addr, $expected, 'construct_address with no params');
}

{
    my $someone_else = Genome::Sys::User->__define__(
        username => join('-', 'fake', Genome::Sys->username, time()),
        email => join('@', 'fake-test-email', Genome::Config::get('email_domain')),
    );
    my $addr = Genome::Utility::Email::construct_address($someone_else->username);
    my $expected = $someone_else->email;
    is($addr, $expected, 'construct_address given a name');
}
