package Genome::Qc::Command::ImportGenotypeVcfFile;

use strict;
use warnings;
use Genome;

class Genome::Qc::Command::ImportGenotypeVcfFile {
    is => 'Genome::Command::DelegatesToResult',
    has => [
        genotype_vcf_file => {
            is => 'Path',
            doc => 'Path to the genotype vcf file to import',
        },
        reference_sequence_build => {
            is => 'Genome::Model::Build::ReferenceSequence',
            id_by => 'reference_sequence_build_id',
            doc => 'Reference sequence build that the genotype vcf is aligned to',
        },
    ],
    has_optional_input => [
        requestor => { via => '__self__', to => 'reference_sequence_build' },
        user => { via => '__self__', to => 'reference_sequence_build' },
        label => { value => 'reference sequence build' },
        reference_sequence_build_id => {
            is => 'Number',
            doc => 'the reference to use by id',
        },
    ],
    doc => 'Import a genotype vcf file for QC. This creates a software result and the genotype vcf file is copied into its allocation.',
};

sub result_class {
    return 'Genome::SoftwareResult::ImportedFile';
}

sub input_hash {
    my $self = shift;
    return (
        source_file_path => $self->genotype_vcf_file,
    );
}

sub post_get_or_create {
    my $self = shift;
    $self->status_message("Software result id: %s", $self->output_result->id);
}

1;
