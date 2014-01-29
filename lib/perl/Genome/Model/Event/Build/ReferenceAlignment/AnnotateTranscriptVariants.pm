package Genome::Model::Event::Build::ReferenceAlignment::AnnotateTranscriptVariants;

use strict;
use warnings;

use Genome;

class Genome::Model::Event::Build::ReferenceAlignment::AnnotateTranscriptVariants{
    is => ['Genome::Model::Event'],
    has => [
        analysis_base_path => {
            doc => "the path at which all analysis output is stored",
            calculate_from => ['build'],
            calculate      => q|
            return $build->variants_directory;
            |,
            is_constant => 1,
        },
        pre_annotation_filtered_snp_file => {
            doc => "",
            calculate_from => ['analysis_base_path'],
            calculate      => q|
            return $analysis_base_path .'/filtered.variants.pre_annotation';
            |,
        },  
        post_annotation_filtered_snp_file => {
            doc => "",
            calculate_from => ['analysis_base_path'],
            calculate      => q|
            return $analysis_base_path .'/filtered.variants.post_annotation';
            |,
        }, 
    ],
};

sub execute {
    my $self = shift;
    my $pre_annot_file  = $self->pre_annotation_filtered_snp_file;
    my $post_annot_file = $self->post_annotation_filtered_snp_file;

    unless ($self->check_for_existence($pre_annot_file)) {
        $self->error_message("Adapted filtered snp file: $pre_annot_file does not exist for annotation");
        return;
    }

    #Allowing empty pre_annotation file
    if (-z $pre_annot_file) {
        $self->warning_message("pre annotation filtered snp file: $pre_annot_file is empty. Skip");
        `touch $post_annot_file`;
        return 1;
    }

    # This has been vetted by the processing profile's _initialize_build method already
    my $annotator_version;
    my $annotator_filter;
    my $accept_reference_IUB_codes;
    {
        my $build = $self->build;
        my $model = $build->model;
        my $pp    = $model->processing_profile;
        $annotator_version = $pp->transcript_variant_annotator_version;
        $annotator_filter  = $pp->transcript_variant_annotator_filter;
        $accept_reference_IUB_codes = $pp->transcript_variant_annotator_accept_reference_IUB_codes;
    }
    

    my %params = (
        variant_file               => $pre_annot_file,
        output_file                => $post_annot_file,
        annotation_filter          => $annotator_filter,
        accept_reference_IUB_codes => $accept_reference_IUB_codes,
        no_headers                 => 1,
        #reference_transcripts => $self->model->annotation_reference_transcripts, 
        use_version                => $annotator_version,
    );
    my $abuild = $self->model->annotation_reference_build;
    $params{build_id} = $abuild->id if $abuild;

    use Data::Dumper;
    $self->debug_message("Annotator params are:\n"
      . Dumper(\%params)
    );

    my $annotator = Genome::Model::Tools::Annotate::TranscriptVariants->create(%params);

    my $rv = $annotator->execute;
    unless ($rv){
        $self->error_message("annotation of adapted filtered snp file failed");
        return;
    }
    return 1;
}

1;
