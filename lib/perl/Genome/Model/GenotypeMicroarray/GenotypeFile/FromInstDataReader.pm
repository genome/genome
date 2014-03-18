package Genome::Model::GenotypeMicroarray::GenotypeFile::FromInstDataReader;

use strict;
use warnings;

use Genome;

class Genome::Model::GenotypeMicroarray::GenotypeFile::FromInstDataReader {
    is => 'Genome::Utility::IO::SeparatedValueReader',
};

sub create {
    my ($class, %params) = @_;

    $params{headers} = []; # holder

    my $self = $class->SUPER::create(%params);
    return if not $self;

    my $resolve_headers = $self->_resolve_headers;
    return if not $resolve_headers;

    return $self;
}

sub _resolve_headers {
    my $self = shift;

    my $header_line;
    do { $header_line = $self->getline; } until not $header_line or $header_line =~ /,/;
    if ( not $header_line ) {
        $self->error_message('Failed to get header line for genotype file!');
        return;
    }

    chomp $header_line;
    my @headers = map { s/\s/_/g; s/_\-\_top$//i; lc } split(',', $header_line);
    $self->headers(\@headers);

    return 1;
}

BEGIN {
    *read = \&Genome::Utility::IO::SeparatedValueReader::next;
}

1;

