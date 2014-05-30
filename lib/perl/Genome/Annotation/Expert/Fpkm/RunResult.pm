package Genome::Annotation::Expert::Fpkm::RunResult;

use strict;
use warnings FATAL => 'all';
use Genome;

class Genome::Annotation::Expert::Fpkm::RunResult {
    is => 'Genome::Annotation::Expert::ResultBase',
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

sub output_filename_base {
    return 'fpkm.vcf';
}

sub output_filename {
    my $self = shift;
    return $self->output_filename_base.'.gz';
}

sub _run {
    my $self = shift;

    Genome::Model::Tools::Vcf::AnnotateWithFpkm->execute(
        vcf_file => $self->input_result_file_path,
        output_file => File::Spec->join($self->temp_staging_directory, $self->output_filename),
        fpkm_file => $self->fpkm_file,
        sample_name => $self->tumor_sample_name,
    );

    return;
}
1;
