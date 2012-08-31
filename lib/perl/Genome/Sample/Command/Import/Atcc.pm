package Genome::Sample::Command::Import::Atcc;

use strict;
use warnings;

use Genome;

class Genome::Sample::Command::Import::Atcc { 
    is => 'Genome::Sample::Command::Import::Base',
    has => [
        name => {
            is => 'Text',
            doc => 'ATCC designation from the provider. This should start with ATCC and have 2 - 3 additional parts separarted by dashes (-).',
        },
        gender => {
            is => 'Text',
            valid_values => [qw/ male female /],
            doc => 'The gender of thepatient.',
        },
        age => {
            is => 'Number',
            doc => 'The age of the patient at time of sample taking.',
        },
        ethnicity => {
            is => 'Text',
            doc => 'The ethnicity of the patient.',
        },
        organ_name => {
            is => 'Text',
            doc => 'The organ from which the sample was taken',
        },
        disease => {
            is => 'Text',
            is_optional => 1,
            doc => 'The name of the disease if present in the sample. If given, the sample common name will be tumor. If no disease is present ther common name will be normal.', 
        },
    ],
};

sub help_brief {
    return 'import ATCC samples';
}

sub execute {
    my $self = shift;

    $self->status_message('Import ATCC Sample...');

    my $name = $self->name;
    my @tokens = split('-', $name);
    if ( $tokens[0] ne 'ATCC' or ( @tokens != 3 and @tokens != 4 ) ) {
        $self->error_message("Failed to validate ATCC name (".$name."). It must start with ATCC and contain 3 -\ 4 parts sparated by dashes.");
        return;
    }

    my %sample_params = (
        name => $name,
        common_name => 'normal',
        extraction_label => $name,
        extraction_type => 'genomic dna',
        cell_type => 'unknown',
        nomenclature => 'ATCC',
        age => $self->age,
        organ_name => $self->organ_name,
    );
    if ( $self->disease ) {
        $sample_params{disease} = $self->disease;
        $sample_params{common_name} = 'tumor';
    }

    my $common_name = join('-', @tokens[1..2]);
    my $individual_name = join('-', @tokens[0..2]);
    my $import = $self->_import(
        taxon => 'human',
        individual => {
            name => $individual_name,
            upn => $individual_name,
            common_name => $common_name,
            nomenclature => 'ATCC',
            gender => $self->gender,
            description => 'ATCC individual: unknown source',
        },
        sample => \%sample_params,
        library => 'extlibs',
    );
    return if not $import;

    $self->status_message('Import...OK');

    return 1;
}

sub _minimum_unique_name_parts {
    return 3;
}

1;

