package Genome::Model::GenotypeMicroarray::Command::Alleles;

use strict;
use warnings;

use Genome;

use Genome::Model::GenotypeMicroarray::GenotypeFile::ReaderForInstData;
require List::Util;

class Genome::Model::GenotypeMicroarray::Command::Alleles {
    is => 'Command::V2',
    has_optional => {
        model => {
            is => 'Genome::Model',
            doc => 'The genotype model to use. This will get the genotype file from the instrument data.',
        },
        build => {
            is => 'Genome::Model::Build',
            doc => 'The genotype build to use. This will get the genotype file from the instrument data.',
        },
        instrument_data => {
            is => 'Genome::InstrumentData',
            doc => 'The genotype instrument data to use. This will get its genotype file.',
        },
    },
    has_optional_transient => {
        alleles => { is => 'Hash', doc => 'Alleles and number of times observed.', },
    },
    has_calculated => {
        total_genotypes => { is => 'Integer', doc => 'Total genotpypes observed.', },
    },
};

sub help_brief { return 'Get the alleles counts from the instrument data genotype file.'; }
sub help_detail { return help_brief(); }

sub total_genotypes {
    my $self = shift;

    my $alleles = $self->alleles;
    return 0 if not $alleles;

    return List::Util::sum( values %$alleles );
}

sub execute {
    my $self = shift;

    my $instdata = $self->_resolve_instrument_data;
    return if not $instdata;

    my $reader = Genome::Model::GenotypeMicroarray::GenotypeFile::ReaderForInstData->create(
        instrument_data => $instdata,
    );
    return if  not $reader;

    my %alleles;
    while ( my $genotype = $reader->read ) {
        $alleles{ $genotype->{allele1} . $genotype->{allele2} }++;
    }
    $self->alleles(\%alleles);

    print "Alleles and observed counts:\n";
    for my $alleles ( sort keys %alleles ) {
        print join(' ', $alleles, $alleles{$alleles})."\n";
    }
    print "Total genotypes: ".$self->total_genotypes."\n";
    
    return 1;
};

sub _resolve_instrument_data {
    my $self = shift;

    my $instdata = $self->instrument_data;
    return $instdata if $instdata;

    if ( $self->build ) {
        return $self->build->instrument_data;
    }

    my $model = $self->model;
    if ( $model ) {
        return $self->build->instrument_data;
    }

    $self->error_message('No model, build, or instrument data given!');
    return;
}

1;

