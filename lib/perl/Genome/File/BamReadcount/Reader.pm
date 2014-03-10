package Genome::File::BamReadcount::Reader;

use Genome::File::TypedStream;
use Genome::File::BamReadcount::Entry;
use Genome;
use Carp qw/confess/;
use strict;
use warnings;

use base qw(Genome::File::TypedStream);

=head1 NAME

Genome::File::BamReadcount::Reader - A class for reading bam-readcount files.

=head1 SYNOPSIS

my $reader = new Genome::File::BamReadcount::Reader("input.rc"); # or input.rc.gz

while (my $entry = $reader->next) {
    # ...
}

$reader->close;

=cut

sub _parse_header {
    return 1;   #bam-readcount files have no header. sad times.
}

sub _next_entry {
    my $self = shift;
    my $line = $self->_getline;
    return unless $line;
    chomp $line;

    my $entry = Genome::File::BamReadcount::Entry->new($line);
    return $entry;
}

1;
