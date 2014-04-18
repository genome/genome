package Genome::Annotation::Joinx::RunResult;

use strict;
use warnings FATAL => 'all';
use Genome;
use File::Spec;

class Genome::Annotation::Joinx::RunResult {
    is => 'Genome::Annotation::Detail::Result',
    has_input => [
        known_variants => {
            is => 'Genome::Model::Build::ImportedVariationList',
            is_many => 1,
        },
    ],
    has_param => [
        info_string => {
            is => 'Text',
        },
        version => {
            is => 'Text',
        },
    ],
};

sub output_filename {
    return 'joinx_vcf_annotate.vcf.gz';
}

sub _run {
    my $self = shift;
    my @known_variants = $self->known_variants;
    if (scalar @known_variants != 1) {
        die "We don't currently support more than one annotation vcf";
    }
    my $annotation_build = $known_variants[0];

    my $input_file  = $self->input_result->output_file_path;
    my $output_file = File::Spec->join($self->temp_staging_directory, $self->output_filename);
    my $info_string = $self->info_string;
    my $info        = $info_string ? 1 : 0;

    my $vcf_annotator = Genome::Model::Tools::Joinx::VcfAnnotate->create(
        input_file      => $input_file,
        annotation_file => $annotation_build->snvs_vcf,
        output_file     => $output_file,
        use_bgzip       => 1,
        info_fields     => $info_string,
        info            => $info,
        use_version     => $self->version,
    );

    unless ($vcf_annotator->execute) {
        die $self->error_message("Failed to execute joinx vcf-annotate");
    }

    return 1;
}


1;
