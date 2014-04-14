package Genome::File::Vcf::Writer;

use Genome;
use Genome::File::Vcf::Header;
use Genome::File::Vcf::Entry;
use Data::Dumper;
use Carp qw/confess/;
use strict;
use warnings;

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
    confess "No vcf header given to vcf writer!" unless $header;
    $fh->print($header->to_string . "\n");

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
        confess "Attempted to write class $type to Vcf file, expected a Vcf entry"
            unless $type eq 'Genome::File::Vcf::Entry';
        $self->{_filehandle}->print($entry->to_string . "\n");
    }
}

sub close {
    my $self = shift;
    $self->{_filehandle}->close;
}

1;
