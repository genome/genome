package Genome::File::Vcf::Reader;

use Data::Dumper;
use Genome;
use Carp qw/confess/;
use strict;
use warnings;

class Genome::File::Vcf::Reader {
    is => "UR::Object",
    has => [
        name => {
            is => "Text",
        },
    ],
    has_transient_optional => [
        filehandle => {
            is => "SCALAR",
        },
        header => {
            is => "Genome::File::Vcf::Header",
        },
        _line_buffer => {
            is => "ARRAY",
        },
    ],
};

sub open {
    my ($self, $path) = @_;
    return $self->fhopen(Genome::Sys->open_file_for_reading($path));
}

sub fhopen {
    my ($self, $fh) = @_;
    $self->filehandle($fh);
    $self->_parse_header();
}

sub _parse_header {
    my $self = shift;
    my @lines;
    my $name = $self->name;
    my $fh = $self->filehandle;

    $self->_line_buffer([]);

    while (my $line = $fh->getline) {
        chomp $line;
        if ($line =~ /^#/) {
            push(@lines, $line);
        } else {
            push(@{$self->_line_buffer}, $line);
            last;
        }
    }
    confess "No vcf header found in file $name" unless @lines;
    my $header = Genome::File::Vcf::Header->create;
    $header->parse(@lines);
    $self->header($header);
}

sub next {
    my $self = shift;
    my $line;
    if (@{$self->_line_buffer}) {
        $line = shift @{$self->_line_buffer};
    } else {
        $line = $self->filehandle->getline;
    }
    chomp $line if $line;
    return unless $line;

    my $entry = Genome::File::Vcf::Entry->create;
    $entry->parse($line);
    return $entry;
}

1;
