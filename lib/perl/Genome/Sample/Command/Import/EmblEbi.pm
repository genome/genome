package Genome::Sample::Command::Import::EmblEbi;

use strict;
use warnings;

use Genome;

class Genome::Sample::Command::Import::EmblEbi { 
    is => 'Genome::Sample::Command::Import::Base',
    has => [
        name => {
            is => 'Text',
            doc => 'EMBL-EBI sample name. This is hyphen separated and includes 3 parts: "EMBL", patient name and sample name. Ex: EMBL-HCT20142-ERS025081.',
        },
        gender => {
            is => 'Text',
            valid_values => [qw/ male female unknown /],
            doc => 'The gender of the patient.',
        },
        age => {
            is => 'Number',
            doc => 'The age of the patient at time of sample taking.',
        },
        ethnicity => {
            is => 'Text',
            doc => 'The ethinicity of the patient.',
        },
        tissue => {
            is => 'Text',
            doc => 'Tissue from where the sample was taken.',
        },
        extraction_type => {
            is => 'Text',
            valid_values => [ 'cdna', 'genomic dna', 'ipr product', 'rna', 'total rna', ],
            doc => 'Type of sample.',
        },
        _patient_name => { is_optional => 1, },
    ],
};

sub help_brief {
    return 'import EMBL-EBI samples';
}

sub execute {
    my $self = shift;

    $self->status_message('Import EMBL-EBI Sample...');

    my $validate_and_set = $self->_validate_name_and_set_individual_name;
    return if not $validate_and_set;

    my $import = $self->_import(
        taxon => 'human',
        individual => {
            name => $self->_patient_name,
            upn => $self->_patient_name,
            nomenclature => 'EMBL-EBI',
            gender => $self->gender,
            ethnicity => $self->ethnicity,
        },
        sample => {
            name => $self->name,
            extraction_label => $self->name,
            extraction_type => $self->extraction_type,
            tissue_desc => $self->tissue,
            tissue_label => $self->tissue,
            cell_type => 'unknown',
            nomenclature => 'EMBL-EBI',
            age => $self->age,
        },
        library => 'extlibs',
    );
    return if not $import;

    $self->status_message('Import...OK');

    return 1;
}
sub _validate_name_and_set_individual_name {
    my $self = shift;

    my $name = $self->name;
    my @tokens = split('-', $name);
    if ( not @tokens == 3 ) {
        $self->error_message("Invalid EMBL-EBI name ($name). It must have 3 parts separated by hyphens.");
        return;
    }

    if ( not $tokens[0] eq 'EMBL' ) {
        $self->error_message("Invalid EMBL-EBI name ($name). It must start with EMBL.");
        return;
    }

    if ( my @invalid_tokens = grep { $_ !~ /^[\w\d]+$/ } @tokens ) {
        $self->error_message("Found invalid characters in EMBL-EBI name. Only letters and numbers are allowed: @invalid_tokens");
        return;
    }
    my $patient_name = join('-', @tokens[0..1]);
    $self->status_message('Patient name: '.$patient_name);
    $self->_patient_name($patient_name);

    return 1
}

1;

