package Genome::File::TypedStream;

use Carp qw/confess/;
use Genome;
use IO::File;

use strict;
use warnings;

sub new {
    my ($class, $filename) = @_;
    my $fh;
    if(Genome::Sys->file_is_gzipped($filename)) {
        $fh = Genome::Sys->open_gzip_file_for_reading($filename);
    }
    else {
        $fh = Genome::Sys->open_file_for_reading($filename);
    }
    return $class->fhopen($fh, $filename);
}

sub fhopen {
    my ($class, $fh, $name) = @_;

    $name |= "Unknown input stream (possibly stdin)";

    my $self = {
        _cached_entry => undef,
        _filehandle => $fh,
        _filters => [],
        _have_putback => 0,
        _last_line => undef,
        header => undef,
        line_number => 0,
        name => $name,
    };

    bless $self, $class;

    $self->_parse_header;

    return $self;
}

sub _getline {
    my $self = shift;

    if ($self->{_have_putback}) {
        $self->{_have_putback} = 0;
        return $self->{_last_line};
    }

    my $line = $self->{_filehandle}->getline;
    $self->{_last_line} = $line;
    if ($line) {
        ++$self->{line_number};
    }
    return $line;
}

sub putback {
    my $self = shift;

    if ($self->{_have_putback}) {
        confess "Attempted to call putback multiple times!";
    }

    $self->{_have_putback} = 1;
}

sub _parse_header {
    my $self = shift;
    confess "_parse_header must be defined!";
}

sub _next_entry {
    my $self = shift;
    confess "_next_entry must be defined!";
}

sub add_filter {
    my ($self, $filter_coderef) = @_;
    push(@{$self->{_filters}}, $filter_coderef);
}

sub header {
    my $self = shift;
    return $self->{header};
}

sub next {
    my $self = shift;
    if ($self->{_cached_entry}) {
        my $cached = $self->{_cached_entry};
        $self->{_cached_entry} = undef;
        return $cached;
    }

    ENTRY: while (my $entry = $self->_next_entry) {
        if (defined $self->{_filters}) {
            for my $filter (@{$self->{_filters}}) {
                next ENTRY unless $filter->($entry);
            }
        }
        return $entry;
    }
    return;
}

sub peek {
    my $self = shift;
    $self->{_cached_entry} = $self->next unless defined $self->{_cached_entry};
    return $self->{_cached_entry};
}

1;
