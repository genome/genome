package Genome::Model::GenotypeMicroarray::GenotypeFile::ReadTsv;

use strict;
use warnings;

use Genome;

class Genome::Model::GenotypeMicroarray::GenotypeFile::ReadTsv { 
    is => 'Genome::Utility::IO::SeparatedValueReader',
    has => {
        separator => { is => 'Text', value => "\t", },
    },
    has_transient => {
        alleles => { is => 'Hash', default_value => {}, },
    },

};

sub read {
    my $self = shift;

    my $genotype = $self->next;
    return if not $genotype;

    $genotype->{id} = $genotype->{snp_name};
    $genotype->{alleles} = $genotype->{allele1}.$genotype->{allele2};

    return $genotype;
}

1;

