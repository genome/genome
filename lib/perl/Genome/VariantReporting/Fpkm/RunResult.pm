package Genome::VariantReporting::Fpkm::RunResult;

use strict;
use warnings FATAL => 'all';
use Genome;

class Genome::VariantReporting::Fpkm::RunResult {
    is => 'Genome::VariantReporting::Framework::Component::Expert::Result',
    has_input => [
        fpkm_file => {
            is => 'String',
        },
        tumor_sample_name => {
            is => 'String',
        },
    ],
    has_param => [
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
        vcf_file => $self->input_result_file_path,
        output_file => $self->fpkm_output_file,
        fpkm_file => $self->fpkm_file,
        sample_name => $self->tumor_sample_name,
    );
    $command->execute;

    return;
}
1;
