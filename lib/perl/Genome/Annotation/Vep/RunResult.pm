package Genome::Annotation::Vep::RunResult;

use strict;
use warnings FATAL => 'all';
use Genome;
use Sys::Hostname;

class Genome::Annotation::Vep::RunResult {
    is => 'Genome::Annotation::Detail::Result',
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
        condel => { is => 'String', },
        terms => {is => 'String', },
        regulatory => {is => 'Boolean',},
        gene => {is => 'Boolean',},
        most_severe => {is => 'Boolean',},
        per_gene => {is => 'Boolean',},
        hgnc => {is => 'Boolean',},
        coding_only => {is => 'Boolean',},
        canonical => {is => 'Boolean',},
        plugins => {is => 'String',
                    is_many => 1},
        plugins_version => {is => 'String',},
        hgvs => {is => 'Boolean',},
    ],
};

sub output_filename_base {
    return 'vep.vcf';
}

sub output_filename {
    my $self = shift;
    return $self->output_filename_base.'.gz';
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

    if ($self->hgvs) {
        $params{fasta} = $self->reference_build->fasta_file;
    }
    my $vep_output_file = File::Spec->join($self->temp_staging_directory, $self->output_filename_base);
    my $final_output_file = File::Spec->join($self->temp_staging_directory, $self->output_filename);
    my $vep_command = Genome::Db::Ensembl::Command::Run::Vep->create(
        input_file => $self->input_result->output_file_path,
        output_file => $vep_output_file,
        ensembl_version => $self->ensembl_version,
        custom => \@custom_annotation_inputs,
        format => "vcf",
        vcf => 1,
        quiet => 0,
        %params,
    );

    unless ($vep_command->execute) {
        die $self->error_message("Failed to execute vep");
    }

    Genome::Sys->gzip_file($vep_output_file, $final_output_file);
    Genome::Model::Tools::Tabix::Index->execute(
        input_file => $final_output_file,
        preset => 'vcf',
    );

    return;
}
