package Genome::Annotation::Vep::Result;

use strict;
use warnings FATAL => 'all';
use Genome;
use Sys::Hostname;

class Genome::Annotation::Vep::Result {
    is => 'Genome::Annotation::Detail::Result',
    has_input => [
        ensembl_version => {
            is => 'String',
        },
        feature_list_ids_and_tags => {
            is => 'String',
            is_many => 1,
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

    my @custom_annotation_inputs;
    for my $feature_list_and_tag ($self->feature_list_ids_and_tags) {
        my ($id, $tag) = split(":", $feature_list_and_tag);
        my $feature_list = Genome::FeatureList->get($id);
        push @custom_annotation_inputs, join("@",
            $feature_list->get_tabix_and_gzipped_bed_file,
            $tag,
            "bed",
            "overlap",
            "0",
        );
    }

    my %params = $self->param_hash;
    delete $params{variant_type};
    delete $params{test_name};

    my $vep_command = Genome::Db::Ensembl::Command::Run::Vep->create(
        input_file => $self->input_vcf_result->output_file_path,
        output_file => File::Spec->join($self->temp_staging_directory, $self->output_filename),
        ensembl_version => $self->ensembl_version,
        custom => \@custom_annotation_inputs,
        %params,
    );

    unless ($vep_command->execute) {
        die $self->error_message("Failed to execute vep");
    }

    return;
}
