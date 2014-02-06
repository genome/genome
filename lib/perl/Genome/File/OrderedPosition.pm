package Genome::File::OrderedPosition;

use Genome;
use IO::File;
use Carp qw/confess/;
use Params::Validate qw/validate_pos :types/;
use IPC::Run qw/run/;

use strict;
use warnings FATAL => "all";

my $HEADER_REGEX = "^chrom\tposition";

sub new {
    my ($class, $filename, $memory_gb) = @_;

    my $header = _read_header($filename);

    my $sorted_filename = Genome::Sys->create_temp_file_path();
    my $sort_memory_gb = $memory_gb / 2;
    run(
        ["grep", "-v", $HEADER_REGEX, $filename], "|",
        ["sort", "-k2,2", "-n", "-S$sort_memory_gb" . "G"], "|",
        ["sort", "-k1,1", "-s", "-S$sort_memory_gb" . "G"],
        ">", "$sorted_filename"
    );
    my $fh = Genome::Sys->open_file_for_reading($sorted_filename);

    my $self = {
        _filehandle => $fh,
        _last_line => undef,
        header => $header,
        line_number => 0,
        _last_chromosome => undef,
        _last_start_position => undef,
        filename => $filename,
        sorted_filename => $sorted_filename,
    };

    bless $self, $class;

    return $self;
}

sub _read_header {
    my $filename = shift;

    my $fh = Genome::Sys->open_file_for_reading($filename);
    my $line = $fh->getline;
    $fh->close;

    if($line =~ m/$HEADER_REGEX/) {
        return $line;
    }
    else {
        return;
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

    while (!$self->reached_position($chr, $start)) {
        return unless $self->getline;
    }

    if ($self->exactly_reached_position($chr, $start)) {
        return $self->{_last_line};
    }
    else {
        return;
    }
}

sub exactly_reached_position {
    my ($self, $chr, $start) = @_;

    if ($self->{_last_chromosome} eq $chr && $self->{_last_start_position} == $start) {
        return 1;
    }
    else {
        return 0;
    }
}

sub reached_position {
    my ($self, $chr, $start) = @_;

    unless (defined($self->{_last_chromosome})) {
        return 0;
    }

    if ($chr gt $self->{_last_chromosome}) {
        return 0;
    }

    if ($chr eq $self->{_last_chromosome} && $start > $self->{_last_start_position}) {
        return 0;
    }

    return 1;
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

sub close {
    my $self = shift;

    $self->{_filehandle}->close;
}

1;
