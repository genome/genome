package Genome::Sample::Command::Import::Metahit;

use strict;
use warnings;

use Genome;

require Carp;
use Data::Dumper 'Dumper';

class Genome::Sample::Command::Import::Metahit { 
    is => 'Genome::Sample::Command::Import::Base',
    has => [
        name => {
            is => 'Text',
            doc => 'MetaHIT sample name.',
        },
        tissue_name => {
            is => 'Text',
            doc => 'Tisse anem from where the sample is from',
        },
        gender => {
            is => 'Text',
            valid_values => [qw/ male female /],
            doc => 'The gender of the individual.',
        },
        age => {
            is => 'Number',
            doc => 'The age of the individual at time of sample taking.',
        },
        bmi => {
            is => 'Text',
            doc => 'The body mass index of the individual at time of sample taking.',
        },
    ],
};

sub help_brief {
    return 'import METAHIT samples';
}

sub execute {
    my $self = shift;

    $self->status_message('Import METAHit Sample...');

    my $individual_name = 'METAHIT-'.$self->name;
    my $sample_name = $individual_name.'-1'; # TODO allow for more than one sample
    my $import = $self->_import(
        taxon => 'human',
        individual => {
            name => $individual_name,
            upn => $individual_name,
            nomenclature => 'METAHIT',
            gender => $self->gender,
            description => 'METAHit individual: unknown source',
        },
        sample => {
            name => $sample_name,
            extraction_label => $sample_name,
            tissue_desc => $self->tissue_name,
            tissue_label => $self->tissue_name,
            extraction_type => 'genomic dna',
            cell_type => 'unknown',
            nomenclature => 'METAHIT',
            age => $self->age,
            body_mass_index => $self->bmi,
        },
        library => 'extlibs',
    );
    return if not $import;

    $self->status_message('Import...OK');

    return 1;
}

1;

