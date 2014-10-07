package Genome::VariantReporting::Expert::Joinx::RunResult;

use strict;
use warnings FATAL => 'all';
use Genome;
use File::Spec;

class Genome::VariantReporting::Expert::Joinx::RunResult {
    is => 'Genome::VariantReporting::Framework::Component::Expert::Result',
    has_input => [
        vcf_lookup => {
            is => 'text',
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
    has_transient_optional => [
        vcf => {
            is => 'Path',
        },
    ],
};


sub output_filename {
    return 'joinx_vcf_annotate.vcf.gz';
}

sub _run {
    my $self = shift;

    my $output_file = File::Spec->join($self->temp_staging_directory, $self->output_filename);
    my $info_string = $self->info_string;
    my $info        = $info_string ? 1 : 0;

    my $vcf_annotator = Genome::Model::Tools::Joinx::VcfAnnotate->create(
        input_file      => $self->input_vcf,
        annotation_file => $self->vcf,
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
