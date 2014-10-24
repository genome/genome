package Genome::File::Breakdancer::Reader;

use Genome::File::Breakdancer::Header;
use Genome::File::Breakdancer::Entry;

use Carp qw/confess/;
use Genome;
use Genome::File::TypedReader;

use strict;
use warnings;

use base qw(Genome::File::TypedReader);

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

    $self->{header} = new Genome::File::Breakdancer::Header(\@lines);
}

sub _next_entry {
    my $self = shift;
    my $line = $self->_getline;
    return unless $line;
    chomp $line;

    my $entry = Genome::File::Breakdancer::Entry->new($self->{header}, $line);
    return $entry;
}

1;
