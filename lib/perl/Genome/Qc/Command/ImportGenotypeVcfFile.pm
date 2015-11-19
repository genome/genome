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
        }
    ],
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

1;
