package Genome::File::BedPe::Writer;

use Carp qw/confess/;
use IO::File;
use Genome;
use Genome::File::BedPe::Header;
use Genome::File::BedPe::Entry;

use strict;
use warnings;

sub new {
    my ($class, $path, $header) = @_;
    return $class->fhopen(new IO::File($path, "w"), $path, $header);
}

sub fhopen {
    my ($class, $fh, $name, $header) = @_;
    $name |= "unknown file path";
    confess "No vcf header given to vcf writer!" unless $header;
    my $str = $header->to_string;
    if ($str) {
        $fh->print("$str\n");
    }

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
        confess "Attempted to write class $type to BedPe file, expected a BedPe entry"
            unless $type eq 'Genome::File::BedPe::Entry';
        $self->{_filehandle}->print($entry->to_string . "\n");
    }
}

1;
