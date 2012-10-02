package Genome::File::Vep::Reader;

use Genome::File::Vep::Header;
use Genome::File::Vep::Entry;
use Data::Dumper;
use Genome;
use Carp qw/confess/;
use strict;
use warnings;

sub new {
    my ($class, $path) = @_;
    return $class->fhopen(Genome::Sys->open_file_for_reading($path), $path);
}

sub fhopen {
    my ($class, $fh, $name) = @_;
    $name |= "unknown vep file";
    my $header_txt = $fh->getline;
    while ($header_txt =~ /^##/) {
        $header_txt = $fh->getline;
    }
    $header_txt =~ s/^#//;

    my $self = {
        name => $name,
        filehandle => $fh,
        header => new Genome::File::Vep::Header($header_txt),
        line_num => 1,
        _cached_entry => undef,
    };

    return bless $self, $class;
}

sub next {
    my $self = shift;
    if ($self->{_cached_entry}) {
        my $rv = $self->{_cached_entry};
        undef $self->{_cached_entry};
        return $rv;
    }
    while (my $line = $self->{filehandle}->getline) {
        ++$self->{line_num};
        chomp $line;
        # There are blank lines in vep files sometimes ._.
        next if !$line || $line =~ /^#/;
        return Genome::File::Vep::Entry->new($line);
    }
}

sub peek {
    my $self = shift;
    $self->{_cached_entry} = $self->next unless $self->{_cached_entry};
    return $self->{_cached_entry};
}

1;
