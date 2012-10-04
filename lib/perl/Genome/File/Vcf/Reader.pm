package Genome::File::Vcf::Reader;

use Genome::File::Vcf::Entry;
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

sub new {
    my ($class, $filename) = @_;
    my $fh;
    if(Genome::Sys->file_is_gzipped($filename)) {
        $fh = Genome::Sys->open_gzip_file_for_reading($filename);
    } else {
        $fh = Genome::Sys->open_file_for_reading($filename);
    }
    return $class->fhopen($fh, $filename);
}

sub fhopen {
    my ($class, $fh, $name) = @_;
    $name |= "unknown vcf file";
    my $self = {
        name => $name,
        filehandle => $fh,
        _header => 0,
        _line_buffer => [],
    };
    bless $self, $class;
    $self->_parse_header();
    return $self;
}

sub _parse_header {
    my $self = shift;
    my @lines;
    my $name = $self->name;
    my $fh = $self->filehandle;

    while (my $line = $fh->getline) {
        chomp $line;
        if ($line =~ /^#/) {
            push(@lines, $line);
        } else {
            push(@{$self->{_line_buffer}}, $line);
            last;
        }
    }
    confess "No vcf header found in file $name" unless @lines;
    $self->{header} = Genome::File::Vcf::Header->create(lines => \@lines);
}

sub next {
    my $self = shift;
    my $line;
    if (@{$self->{_line_buffer}}) {
        $line = shift @{$self->{_line_buffer}};
    } else {
        $line = $self->{filehandle}->getline;
    }
    chomp $line if $line;
    return unless $line;

    my $entry = Genome::File::Vcf::Entry->new($self->{header}, $line);
    return $entry;
}

sub header {
    my $self = shift;
    return $self->{header};
}

1;
