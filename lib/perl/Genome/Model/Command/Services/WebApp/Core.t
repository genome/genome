#!/usr/bin/env genome-perl

use strict;
use warnings;

use Test::More tests => 3;
use DBI;
use above 'Genome'; #Otherwise we're testing the deployed code instead when running the test by hand

our @dbi_classes = ();
## I want to trap any call to DBI->connect

sub my_connect {
    my $i = 1;
    while (1) {
        my @caller = caller($i++);
        last unless @caller;

#        print $caller[0], $caller[3], "\n";

        if ($caller[0] =~ /^Genome/) {
            push @dbi_classes, $caller[0];
            last;
        }
    }

    shift->real_connect(@_);
}

*DBI::real_connect = *DBI::connect;
*DBI::connect = *main::my_connect;

0 && &DBI::real_connect; # prevents warning

use_ok('Genome::Model::Command::Services::WebApp::Core');

foreach my $c (@Genome::Model::Command::Services::WebApp::Core::error_classes) {
    diag('Error loading: ' . $c);
}
ok(scalar @Genome::Model::Command::Services::WebApp::Core::error_classes == 0, 'no errors loading classes');

foreach my $c (@dbi_classes) {
    diag('DBI->connect while loading: ' . $c);
}
ok(scalar @dbi_classes == 0, 'no classes caused DBI connect during loading');
