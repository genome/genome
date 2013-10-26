package Genome::File::Vep::Reader;

use Genome::File::TypedStream;
use Genome::File::Vep::Header;
use Genome::File::Vep::Entry;
use Data::Dumper;
use Genome;
use Carp qw/confess/;
use strict;
use warnings;

use base qw(Genome::File::TypedStream);

sub _parse_header {
    my $self = shift;

    my $header_txt = $self->_getline;
    my @lines;
    while ($header_txt =~ /^##/) {
        push(@lines, $header_txt);
        $header_txt = $self->_getline;
    }
    push(@lines, $header_txt);

    $self->{header} = new Genome::File::Vep::Header(\@lines);
}

sub _next_entry {
    my $self = shift;
    while (my $line = $self->_getline) {
        chomp $line;
        # There are blank lines in vep files sometimes ._.
        next if !$line || $line =~ /^#/;
        return Genome::File::Vep::Entry->new($line);
    }
}

1;
