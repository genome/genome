package Genome::Model::GenotypeMicroarray::GenotypeFile::FromBuildOriginalTsvReader;

use strict;
use warnings;

use parent 'Genome::Model::GenotypeMicroarray::GenotypeFile::FromBaseReader';

sub source_type {
    return 'build';
}

sub get_genotype_file {
    my $self = shift;
    return $self->{build}->original_genotype_file_path;
}

sub _resolve_headers {
    my $self = shift;

    my $header_line = $self->{genotype_fh}->getline;
    if ( not $header_line ) {
        $self->error_message('Failed to get header line for genotype file!');
        return;
    }

    chomp $header_line;
    $self->{headers} = [ split($self->{delimiter}, $header_line) ];

    return 1;
}

1;

