package Genome::File::TypedWriter;

use strict;
use warnings;

use Genome;
use Carp qw/confess/;

sub new {
    my ($class, $path, $header) = @_;
    my $fh;
    if ($path =~ m/\.gz$/) {
        $fh = Genome::Sys->open_gzip_file_for_writing($path);
    } else {
        $fh = Genome::Sys->open_file_for_writing($path);
    }
    return $class->fhopen($fh, $path, $header);
}

sub fhopen {
    my ($class, $fh, $name, $header) = @_;
    confess "No header given to writer object ($class)!" unless $header;
    my $str = $header->to_string;
    $fh->print("$str\n") if $str;

    my $self = {
        name => $name,
        _filehandle => $fh,
    };

    return bless $self, $class;
}

sub write {
    my ($self, @entries) = @_;
    for my $entry (@entries) {
        my $type = ref $entry;
        $self->{_filehandle}->print($entry->to_string . "\n");
    }
}

sub close {
    my $self = shift;
    $self->{_filehandle}->close;
}

1;

