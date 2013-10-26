package Genome::File::BedPe::Reader;

use Genome::File::BedPe::Header;
use Genome::File::BedPe::Entry;

use Carp qw/confess/;
use Genome;
use Genome::File::TypedStream;

use strict;
use warnings;

use base qw(Genome::File::TypedStream);

sub _parse_header {
    my $self = shift;
    my @lines;
    my $name = $self->{name};
    my $fh = $self->{_filehandle};

    while (my $line = $self->_getline) {
        chomp $line;
        if ($line =~ /^#/) {
            push(@lines, $line);
        } else {
            $self->putback;
            last;
        }
    }

    $self->{header} = new Genome::File::BedPe::Header(\@lines);
}

sub _next_entry {
    my $self = shift;
    my $line = $self->_getline;
    return unless $line;
    chomp $line;

    my $entry = Genome::File::BedPe::Entry->new($self->{header}, $line);
    return $entry;
}

1;
