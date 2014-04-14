package Genome::Annotation::Vep::Result;

use strict;
use warnings FATAL => 'all';
use Genome;
use Sys::Hostname;

class Genome::Annotation::Vep::Result {
    is => 'Genome::Annotation::Detail::Result',
    has_input => [
        ensembl_annotation_build_id => {
            is => 'String',
        },
        target_region_set => {
            is => 'Genome::FeatureList',
        },
        segmental_duplications_list => {
            is => 'Genome::FeatureList',
        },
        input_vcf_result => {
            is => 'Genome::SoftwareResult',
        },
    ],
    has_param => [
        variant_type => { is => 'Text', },
        format => { is => 'String', },
        polyphen => { is => 'String', },
        sift => { is => 'String', },
        condel => { is => 'String', },
        quiet => { is => 'String', },
    ],
};

sub output_filename {
    return 'vep.vcf';
}

sub output_file_path {
    my $self = shift;
    return File::Spec->join($self->output_dir, $self->output_filename);
}

sub _run {
    my $self = shift;

    my $roi_input = join("@", $self->target_region_set->get_tabix_and_gzipped_bed_file,
        "ROI",
        "bed",
        "overlap",
        "0",
    );
    my $segdup_input = join("@", $self->segmental_duplications_list->get_tabix_and_gzipped_bed_file,
        "SEGDUP",
        "bed",
        "overlap",
        "0",
    );
    my @custom_annotation_inputs = ($roi_input, $segdup_input);

    my %params = $self->param_hash;
    delete $params{variant_type};
    delete $params{test_name};

    my $vep_command = Genome::Db::Ensembl::Command::Vep->create(
        input_file => $self->input_vcf_result->output_file_path,
        output_file => File::Spec->join($self->temp_staging_directory, $self->output_filename),
        ensembl_annotation_build_id => $self->ensembl_annotation_build_id,
        custom => \@custom_annotation_inputs,
        %params,
    );

    unless ($vep_command->execute) {
        die $self->error_message("Failed to execute vep");
    }

    return;
}
