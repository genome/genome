package Genome::File::OrderedPosition;

use Genome;
use IO::File;
use Carp qw/confess/;
use Params::Validate qw/validate_pos :types/;
use IPC::Run qw/run/;

use strict;
use warnings;

sub new {
    my ($class, $filename, $memory) = @_;
    my $temp_file = Genome::Sys->create_temp_file_path();
    run(["sort", "-k1,1", "-k2,2", "-n", "-S$memory", $filename,], ">", "$temp_file");
    my $fh = Genome::Sys->open_file_for_reading($temp_file);
    return $class->fhopen($fh, $filename, $temp_file);
}

sub fhopen {
    my ($class, $fh, $filename, $sorted_filename) = @_;

    my $name |= "Unknown input stream (possibly stdin)";

    my $self = {
        _filehandle => $fh,
        _last_line => undef,
        header => undef,
        line_number => 0,
        _last_chromosome => undef,
        _last_start_position => undef,
        filename => $filename,
        sorted_filename => $sorted_filename,
    };

    bless $self, $class;

    $self->_read_header;

    return $self;
}

sub _read_header {
    my $self = shift;

    my $line = $self->{_filehandle}->getline;

    if($line =~ m/^chrom\tposition/) {
        $self->{header} = $line;
        $self->getline();
    }
    else {
        $self->_parse_line($line);
    }
}

sub getline {
    my $self = shift;

    my $line = $self->{_filehandle}->getline;
    if (defined($line)) {
        return $self->_parse_line($line);
    }
    else {
        return;
    }
}

sub getline_for_position {
    my ($self, $chr, $start) = validate_pos(
        @_, {type => OBJECT}, {type => SCALAR}, {type => SCALAR}
    );

    while (
        $self->{_last_chromosome} <= $chr &&
        $self->{_last_start_position} < $start
    ) {
        return unless $self->getline;
    }

    if (
        $self->{_last_chromosome} == $chr &&
        $self->{_last_start_position} == $start
    ) {
        return $self->{_last_line};
    }
    else {
        return;
    }
}

sub _parse_line {
    my ( $self, $line ) = @_;

    if (defined($line)) {
        my ($chr, $start) = split(/\t/, $line);
        $self->{_last_chromosome} = $chr;
        $self->{_last_start_position} = $start;
        $self->{line_number}++;
        $self->{_last_line} = $line;
        return $line;
    }
    else {
        return;
    }
}

1;
