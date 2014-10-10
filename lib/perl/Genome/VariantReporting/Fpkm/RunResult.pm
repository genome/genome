package Genome::VariantReporting::Fpkm::RunResult;

use strict;
use warnings FATAL => 'all';
use Genome;

class Genome::VariantReporting::Fpkm::RunResult {
    is => 'Genome::VariantReporting::Framework::Component::Expert::Result',
    has_input => [
        fpkm_file_lookup => {
            is => 'Text',
        },
        sample_name => {
            is => 'String',
        },
    ],
    has_param => [
    ],
    has_transient_optional => [
        fpkm_file => {
            is => 'Path',
        },
    ],
};

sub output_filename {
    return 'fpkm.vcf.gz';
}

sub fpkm_output_file {
    my $self = shift;
    return File::Spec->join($self->temp_staging_directory, $self->output_filename);
}

sub _run {
    my $self = shift;

    my $command = Genome::Model::Tools::Vcf::AnnotateWithFpkm->create(
        vcf_file => $self->input_vcf,
        output_file => $self->fpkm_output_file,
        fpkm_file => $self->fpkm_file,
        sample_name => $self->sample_name,
    );
    $command->execute;

    return;
}
1;
