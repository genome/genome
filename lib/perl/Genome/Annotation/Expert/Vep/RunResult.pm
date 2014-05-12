package Genome::Annotation::Expert::Vep::RunResult;

use strict;
use warnings FATAL => 'all';
use Genome;
use Sys::Hostname;

class Genome::Annotation::Expert::Vep::RunResult {
    is => 'Genome::Annotation::Expert::ResultBase',
    has_input => [
        ensembl_version => {
            is => 'String',
        },
        feature_list_ids_and_tags => {
            is => 'String',
            is_many => 1,
        },
        species => {
            is => 'String',
        },
        reference_build => {
            is => 'Genome::Model::Build::ReferenceSequence',
        },
    ],
    has_param => [
        polyphen => { is => 'String', },
        sift => { is => 'String', },
        terms => {is => 'String', },
        regulatory => {is => 'Boolean',},
        canonical => {is => 'Boolean',},
        plugins => {is => 'String',
                    is_many => 1},
        plugins_version => {is => 'String',},
    ],
};

my $BUFFER_SIZE = '5000';

sub output_filename_base {
    return 'vep.vcf';
}

sub output_filename {
    my $self = shift;
    return $self->output_filename_base.'.gz';
}

sub _run {
    my $self = shift;


    my $vep_command = Genome::Db::Ensembl::Command::Run::Vep->create(
        output_file => $self->vep_output_file,
        $self->vep_params,
    );

    unless ($vep_command->execute) {
        die $self->error_message("Failed to execute vep");
    }

    Genome::Sys->gzip_file($self->vep_output_file, $self->final_output_file);
    unlink $self->vep_output_file;

    return;
}

sub vep_output_file {
    my $self = shift;
    return File::Spec->join($self->temp_staging_directory, $self->output_filename_base);
}

sub final_output_file {
    my $self = shift;
    return File::Spec->join($self->temp_staging_directory, $self->output_filename);
}

sub custom_annotation_inputs {
    my $self = shift;

    my $result = [];
    for my $feature_list_and_tag ($self->feature_list_ids_and_tags) {
        my ($id, $tag) = split(":", $feature_list_and_tag);
        my $feature_list = Genome::FeatureList->get($id);
        push @$result, join("@",
            $feature_list->get_tabix_and_gzipped_bed_file,
            $tag,
            "bed",
            "overlap",
            "0",
        );
    }
    return $result;
}

sub vep_params {
    my $self = shift;

    my %params = (
        $self->param_hash,
        input_file => $self->input_result_file_path,
        fasta => $self->reference_build->fasta_file,
        ensembl_version => $self->ensembl_version,
        custom => $self->custom_annotation_inputs,
        format => "vcf",
        vcf => 1,
        quiet => 0,
        hgvs => 1,
        pick => 1,
        buffer_size => $BUFFER_SIZE,
    );
    delete $params{variant_type};
    delete $params{test_name};

    return %params;
}
