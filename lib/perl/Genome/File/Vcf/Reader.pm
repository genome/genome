package Genome::File::Vcf::Reader;

use Genome::File::Vcf::Entry;
use Genome;
use Carp qw/confess/;
use strict;
use warnings;

=head1 NAME

Genome::File::Vcf::Reader - A class for reading vcf files.

=head1 SYNOPSIS

my $reader = new Genome::File::Vcf::Reader("input.vcf"); # or input.vcf.gz
my $header = $reader->header;

while (my $entry = $reader->next) {
    # ...
}

$reader->close;

=cut

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
        _filehandle => $fh,
        _header => undef,
        _line_buffer => [],
        _filters => [],
    };
    bless $self, $class;
    $self->_parse_header();
    return $self;
}

sub close {
    my $self = shift;
    $self->{_filehandle}->close();
}

sub add_filter {
    my ($self, $filter_coderef) = @_;
    push(@{$self->{_filters}}, $filter_coderef);
}

sub _parse_header {
    my $self = shift;
    my @lines;
    my $name = $self->{name};
    my $fh = $self->{_filehandle};

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
    $self->{_header} = Genome::File::Vcf::Header->create(lines => \@lines);
}

sub _next_entry {
    my $self = shift;
    my $line;
    if (@{$self->{_line_buffer}}) {
        $line = shift @{$self->{_line_buffer}};
    } else {
        $line = $self->{_filehandle}->getline;
    }
    chomp $line if $line;
    return unless $line;

    my $entry = Genome::File::Vcf::Entry->new($self->{_header}, $line);
    return $entry;
}

sub next {
    my $self = shift;
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

sub header {
    my $self = shift;
    return $self->{_header};
}

1;
