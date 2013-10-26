package Genome::File::Xsv;
use strict;
use warnings;
use Genome;

class Genome::File::Xsv {
    is => 'Genome::File::Base',
    is_abstract => 1,
    doc => "a file with one row per record, one column per field, with separated values",
};

sub create_reader {
    my $self = shift;
    my @args = @_;
    my $reader = Genome::Utility::IO::SeparatedValueReader->create(separator => $self->separator, input => $self->id, @args);
    die "failed to create writer! " .  Genome::Utility::IO::SeparatedValueReader->error_message() unless $reader;
    return $reader;
}

sub create_writer {
    my $self = shift;
    my @args = @_;
    my $writer = Genome::Utility::IO::SeparatedValueWriter->create(separator => $self->separator, output => $self->id, @args);
    die "failed to create writer! " .  Genome::Utility::IO::SeparatedValueWriter->error_message() unless $writer;
    return $writer;
}

1;

