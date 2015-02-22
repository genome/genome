package Genome::Model::Tools::Sx::PhredReaderBase;

use strict;
use warnings;

use Genome;

class Genome::Model::Tools::Sx::PhredReaderBase {
    is => 'Genome::Model::Tools::Sx::SeqReader',
    is_abstract => 1,
};

sub _parse_io {
    my $self = shift;

    local $/ = "\n>";

    my $entry = $self->{_file}->getline;
    return unless defined $entry;
    chomp $entry;

    if ( $entry eq '>' )  { # very first one
        $entry = $self->{_file}->getline;
        return unless $entry;
        chomp $entry;
    }

    my ($header, $data) = split(/\n/, $entry, 2);
    defined $data && $data =~ s/>//g;

    my ($id, $desc) = split(/\s+/, $header, 2);
    if ( not defined $id or $id eq '' ) {
        Carp::confess("Cannot get id from header ($header) for entry:\n$entry");
    }
    $id =~ s/>//;

    if ( not defined $data ) {
        Carp::confess("No data found for $id entry:\n$entry");
    }

    return ($id, $desc, $data);
}

sub close {
    my $self = shift;

    $self->{_file}->close if $self->{_file};

    return 1;
}

1;

