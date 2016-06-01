package Genome::Model::ReferenceAlignment::Command::AnnotateVariants;

use strict;
use warnings;

use File::Spec;

use Genome;

class Genome::Model::ReferenceAlignment::Command::AnnotateVariants {
    is => 'Genome::Model::ReferenceAlignment::Command::PipelineBase',
    has_param => [
        lsf_queue => {
            default => Genome::Config::get('lsf_queue_build_worker_alt'),
        },
        lsf_resource => {
            default => Genome::Config::get('lsf_resource_annotate_variants'),
        },
    ],
};

sub shortcut {
    my $self = shift;

    return $self->_should_skip_run;
}

sub execute {
    my $self = shift;
    my $build = $self->build;

    return 1 if $self->_should_skip_run;

    my $pp = $build->processing_profile;

    my $pre_annotation_file = $self->_prepare_pre_annotation_file;
    my $post_annotation_file = File::Spec->join($build->variants_directory, 'filtered.variants.post_annotation');

    #Allowing empty pre_annotation file
    if (-z $pre_annotation_file) {
        $self->warning_message("pre annotation filtered snp file: $pre_annotation_file is empty. Skip");
        `touch $post_annotation_file`;
        return 1;
    }

    my $annotator_version = $pp->transcript_variant_annotator_version;
    my $annotator_filter = $pp->transcript_variant_annotator_filter;
    my $accept_reference_IUB_codes = $pp->transcript_variant_annotator_accept_reference_IUB_codes;

    my %params = (
        variant_file               => $pre_annotation_file,
        output_file                => $post_annotation_file,
        annotation_filter          => $annotator_filter,
        accept_reference_IUB_codes => $accept_reference_IUB_codes,
        no_headers                 => 1,
        #reference_transcripts => $self->model->annotation_reference_transcripts, 
        use_version                => $annotator_version,
    );
    my $abuild = $self->build->annotation_reference_build;
    $params{build_id} = $abuild->id if $abuild;

    use Data::Dumper;
    $self->debug_message("Annotator params are:\n"
      . Dumper(\%params)
    );

    my $annotator = Genome::Model::Tools::Annotate::TranscriptVariants->create(%params);

    my $rv = $annotator->execute;
    unless ($rv){
        $self->fatal_message("annotation failed");
    }

    return 1;
}

sub _prepare_pre_annotation_file {
    my $self = shift;
    my $build = $self->build;
    my $pp = $build->processing_profile;

    my $annotator_version = $pp->transcript_variant_annotator_version;
    my $adaptor_version = "Genome::Model::Tools::Annotate::TranscriptVariants::Version" . $annotator_version . "::BedToAnnotation";

    my %adaptor_params;
    $adaptor_params{snv_file} = $build->get_variant_bed_file('snvs', 'hq') if defined $pp->snv_detection_strategy;
    $adaptor_params{indel_file} = $build->get_variant_bed_file('indels', 'hq') if defined $pp->indel_detection_strategy;

    my $pre_annotation_file = File::Spec->join($build->variants_directory, 'filtered.variants.pre_annotation');

    my $adaptor = $adaptor_version->create(
        output => $pre_annotation_file,
        %adaptor_params,
    );
    unless ($adaptor) {
        $self->fatal_message("Could not create annotation adaptor object!");
    }

    my $rv = $adaptor->execute;
    unless ($rv == 1) {
        $self->fatal_message("Problem executing annotation adaptor!");
    }

    unless( Genome::Sys->check_for_path_existence($pre_annotation_file) ){
        $self->error_message("filtered variant output file from find variations step doesn't exist");
        return;
    }

    return $pre_annotation_file;
}

sub _should_skip_run {
    my $self = shift;
    my $build = $self->build;
    my $pp = $build->processing_profile;

    if(($pp->snv_detection_strategy || $pp->indel_detection_strategy) and
            $build->annotation_reference_build) {
        return;
    }
    else {
        $self->debug_message('Decided to skip annotation step.');
        return 1;
    }

}

1;
