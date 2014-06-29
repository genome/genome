package Genome::VariantReporting::Expert::Dbsnp::RunResult;

use strict;
use warnings FATAL => 'all';
use Genome;
use File::Spec;

class Genome::VariantReporting::Expert::Dbsnp::RunResult {
    is => 'Genome::VariantReporting::Component::Expert::Result',
    has_input => [
        known_variants => {
            is => 'Genome::Model::Build::ImportedVariationList',
        },
    ],
    has_param => [
        info_string => {
            is => 'Text',
        },
        joinx_version => {
            is => 'Text',
        },
    ],
};

sub output_filename {
    return 'joinx_vcf_annotate.vcf.gz';
}

sub _run {
    my $self = shift;

    my $input_file  = $self->input_result_file_path;
    my $output_file = File::Spec->join($self->temp_staging_directory, $self->output_filename);
    my $info_string = $self->info_string;
    my $info        = $info_string ? 1 : 0;

    my $vcf_annotator = Genome::Model::Tools::Joinx::VcfAnnotate->create(
        input_file      => $input_file,
        annotation_file => $self->known_variants->snvs_vcf,
        output_file     => $output_file,
        use_bgzip       => 1,
        info_fields     => $info_string,
        info            => $info,
        use_version     => $self->joinx_version,
    );

    unless ($vcf_annotator->execute) {
        die $self->error_message("Failed to execute joinx vcf-annotate");
    }

    return 1;
}


1;
