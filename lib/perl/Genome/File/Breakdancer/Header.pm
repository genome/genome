package Genome::File::Breakdancer::Header;

use strict;
use warnings;

sub new {
    my ($class, $lines) = @_;

    $lines = [] unless defined $lines;

    my $self = {
        lines => $lines,
    };

    bless $self, $class;
    return $self;
}

sub to_string {
    my $self = shift;

    return join("\n", @{$self->{lines}});
}

1;
