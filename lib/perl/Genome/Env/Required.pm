package Genome::Env::Required;

use strict;
use warnings;

use Carp;

sub import {
    my $package = caller;
    return unless $package =~ m/^Genome::Env::/;

    my($env_name) = $package =~ m/::(\w+)$/;

    unless ($ENV{$env_name}) {
        Carp::croak("Environment variable $env_name must be set in your environment or by a site configuration module");
    }
}

1;

=pod

=head1 NAME

Genome::Env::Required

=head1 DESCRIPTION

Make your enviroment variable's module inherit from this one, and the 
enviroment variable becomes required.  If not set, it will throw an exception.

=cut
