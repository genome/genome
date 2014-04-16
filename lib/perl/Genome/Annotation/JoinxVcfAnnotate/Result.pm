package Genome::Annotation::JoinxVcfAnnotate::Result;

use strict;
use warnings FATAL => 'all';
use Genome;
use File::Spec;

class Genome::Annotation::JoinxVcfAnnotate::Result {
    is => 'Genome::Annotation::Detail::Result',
    has_input => [
        annotation_builds => {
            is => 'Genome::Model::Build::ImportedVariationList',
            is_many => 1,
        },
        input_vcf_result => {
            is => 'Genome::SoftwareResult',
        },
    ],
    has_param => [
        variant_type => { 
            is => 'Text', 
        },
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

sub output_file_path {
    my $self = shift;
    return File::Spec->join($self->output_dir, $self->output_filename);
}

sub _run {
    my $self = shift;
    my @annotation_builds = $self->annotation_builds;
    if (scalar @annotation_builds != 1) {
        die "We don't currently support more than one annotation vcf";
    }
    my $annotation_build = $annotation_builds[0];

    my $input_file  = $self->input_vcf_result->output_file_path;
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
        use_version     => $self->joinx_version,
    );

    unless ($vcf_annotator->execute) {
        die $self->error_message("Failed to execute joinx vcf-annotate");
    }

    return 1;
}

