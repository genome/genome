use strict;
use warnings;

package Genome::Env::Required;
use base 'Genome::Env';

sub import {
    my $class = shift;

    return if ($class eq __PACKAGE__);

    my $NAME = $class->NAME();
    unless ($ENV{$NAME}) {
        warn("Environment variable $NAME must be configured to use Genome");
        exit 255;
    }

    return 1;
}

1;

__END__

=pod

=head1 NAME

Genome::Env::Required

=head1 DESCRIPTION

Make your enviroment variable's module inherit from this one, and the
enviroment variable becomes required.  If not set, it will throw an exception.

=cut
