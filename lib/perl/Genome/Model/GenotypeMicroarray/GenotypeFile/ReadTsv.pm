package Genome::Model::GenotypeMicroarray::GenotypeFile::ReadTsv;

use strict;
use warnings;

use Genome;

class Genome::Model::GenotypeMicroarray::GenotypeFile::ReadTsv { 
    is => 'Genome::Utility::IO::SeparatedValueReader',
    has => {
        separator => { is => 'Text', value => "\t", },
    },
};

sub read {
    my $self = shift;

    my $genotype = $self->next;
    return if not $genotype;

    $genotype->{alleles} = $genotype->{allele1}.$genotype->{allele2};

    return $genotype;
}

1;

