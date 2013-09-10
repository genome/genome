package Genome::File::Vcf::Reader;

use Genome::File::TypedStream;
use Genome::File::Vcf::Entry;
use Genome;
use Carp qw/confess/;
use strict;
use warnings;

use base qw(Genome::File::TypedStream);

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

sub _parse_header {
    my $self = shift;
    my @lines;
    my $name = $self->{name};
    my $fh = $self->{_filehandle};

    while (my $line = $self->_getline) {
        chomp $line;
        if ($line =~ /^##/) {
            push(@lines, $line);
        } elsif ($line =~ /^#/) {
            push(@lines, $line);
            last;
        }
        else {
            confess "Invalid vcf header in file $name at line $self->{line_number}";
        }
    }
    confess "No vcf header found in file $name" unless @lines;
    $self->{header} = Genome::File::Vcf::Header->create(lines => \@lines);
}

sub _next_entry {
    my $self = shift;
    my $line = $self->_getline;
    return unless $line;
    chomp $line;

    my $entry = Genome::File::Vcf::Entry->new($self->{header}, $line);
    return $entry;
}

1;
