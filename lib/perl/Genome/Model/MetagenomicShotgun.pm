package Genome::Model::MetagenomicShotgun;

use strict;
use warnings;

use Genome;

class Genome::Model::MetagenomicShotgun {
    is => 'Genome::ModelDeprecated',
    has_param => [
        filter_contaminant_fragments => {
            is => 'Boolean',
            doc => 'when set, reads with mate mapping to contamination reference are considered contaminated, and not passed on to subsequent alignments',
            default => 0,
        },
        filter_duplicates => {
            is => 'Boolean',
            doc =>  'when set, duplicate reads are filtered out when extracting unaligned reads from contamination screen alignment',
        },
        contamination_screen_pp_id => {
            is => 'Text',
            doc => 'processing profile id to use for contamination screen',
            is_optional=> 1,
        },
        metagenomic_nucleotide_pp_id => {
            is => 'Text',
            doc => 'processing profile id to use for metagenomic alignment',
        },
        metagenomic_protein_pp_id => {
            is => 'Text',
            doc => 'processing profile id to use for realignment of unaligned reads from first metagenomic alignment',
            is_optional => 1,
        },
        viral_nucleotide_pp_id => {
            is => 'Text',
            doc => 'processing profile id to use for first viral verification alignment',
            is_optional => 1,
        },
        viral_protein_pp_id => {
            is => 'Text',
            doc => 'processing profile id to use for first viral verification alignment',
            is_optional => 1,
        },
    ],
    has => [
        contamination_screen_reference => {
            is => 'Genome::Model::Build::ImportedReferenceSequence',
            is_mutable => 1,
            is_optional => 1,
            via => 'inputs',
            to => 'value',
            where => [name => 'contamination_screen_reference', value_class_name => 'Genome::Model::Build::ImportedReferenceSequence'],
        },
        metagenomic_protein_reference => {
            is => 'Genome::Model::Build::ImportedReferenceSequence',
            is_mutable => 1,
            via => 'inputs',
            to => 'value',
            where => [name => 'metagenomic_protein_reference ', value_class_name => 'Genome::Model::Build::ImportedReferenceSequence'],
        },
        viral_nucleotide_reference => {
            is => 'Genome::Model::Build::ImportedReferenceSequence',
            is_mutable => 1,
            is_optional => 1,
            via => 'inputs',
            to => 'value',
            where => [name => 'viral_nucleotide_reference ', value_class_name => 'Genome::Model::Build::ImportedReferenceSequence'],
        },
        viral_protein_reference => {
            is => 'Genome::Model::Build::ImportedReferenceSequence',
            is_mutable => 1,
            is_optional => 1,
            via => 'inputs',
            to => 'value',
            where => [name => 'viral_protein_reference ', value_class_name => 'Genome::Model::Build::ImportedReferenceSequence'],
        },
        metagenomic_nucleotide_reference => {
            is => 'Genome::Model::Build::ImportedReferenceSequence',
            is_mutable => 1,
            via => 'inputs',
            to => 'value',
            where => [name => 'metagenomic_alignment_reference', value_class_name => 'Genome::Model::Build::ImportedReferenceSequence'],
        },
        contamination_screen_model => {
            is => 'Genome::Model::ReferenceAlignment',
            is_optional => 1,
            via => 'from_model_links',
            to => 'from_model',
            where => [role => 'contamination_screen_model'],
        },
        metagenomic_nucleotide_model => {
            is => 'Genome::Model::ReferenceAlignment',
            via => 'from_model_links', 
            to => 'from_model',
            where => [role => 'metagenomic_nucleotide_model'],
        },
        metagenomic_protein_model => {
            is => 'Genome::Model::ReferenceAlignment',
            is_optional => 1,
            via => 'from_model_links',
            to => 'from_model',
            where => [role => 'metagenomic_protein_model'],
        },
        viral_nucleotide_model => {
            is => 'Genome::Model::ReferenceAlignment',
            is_optional => 1, 
            via => 'from_model_links',
            to => 'from_model',
            where => [role => 'viral_nucleotide_model'],
        },
        viral_protein_model => {
            is => 'Genome::Model::ReferenceAlignment',
            is_optional => 1, 
            via => 'from_model_links',
            to => 'from_model',
            where => [role => 'viral_protein_model'],
        },
    ],
};

sub sub_model_labels {
    return (qw/ contamination_screen metagenomic_nucleotide metagenomic_protein viral_nucleotide viral_protein /);
}

sub build_subclass_name {
    return 'metagenomic-composition-shotgun';
}

sub delete {
    my $self = shift;
    for my $sub_model ($self->from_models) {
        $sub_model->delete;
    }
    return $self->SUPER::delete(@_);
}

sub create {
    my ($class, %params) = @_;

    my $self = $class->SUPER::create(%params);
    return unless $self;

    my $processing_profile = $self->processing_profile;
    for ( $self->sub_model_labels ){
        my $pp_method = $_."_pp_id";
        if($self->processing_profile->$pp_method) {
            my $model = $self->_create_model_for_type($_);
            unless ($model) {
                $self->error_message("Error creating $_ model!");
                $self->delete;
                return;
            }
        }
    }

    return $self;
}

sub sequencing_platform{
    return 'solexa';
}

# TODO make this get or create
sub _create_model_for_type {
    my $self = shift;
    my $type = shift;

    #CREATE UNDERLYING REFERENCE ALIGNMENT MODELS
    my $pp_accessor = $type."_pp_id";
    my $reference_accessor = $type."_reference";
    my %model_params = (
        processing_profile_id => $self->processing_profile->$pp_accessor,
        subject_name => $self->subject_name,
        name => $self->name.".$type model",
        reference_sequence_build=>$self->$reference_accessor,
    );
    my $model = Genome::Model::ReferenceAlignment->create( %model_params );

    unless ($model){
        die $self->error_message("Couldn't create contamination screen model with params ".join(", ", map {$_ ."=>". $model_params{$_}} keys %model_params) );
    }

    $self->add_from_model(from_model=> $model, role=>$type.'_model');
    $self->status_message("Created $type model ".$model->__display_name__);

    return $model;
}

sub _resolve_resource_requirements_for_build {
    my ($self, $build) = @_;
    my @instrument_data = $build->instrument_data;
    my $gtmp = 30 + 5 * (1 + scalar(@instrument_data));
    return "-R 'rusage[gtmp=$gtmp:mem=16000]' -M 16000000";
}

sub _execute_build {
    my ($self, $build) = @_;

    my $model = $build->model;
    $self->status_message('Build '.$model->__display_name__);
    my $contamination_screen_model = $model->contamination_screen_model;
    $self->status_message("Got contamination_screen_model ".$contamination_screen_model->__display_name__) if $contamination_screen_model;
    my $metagenomic_nucleotide_model = $model->metagenomic_nucleotide_model;
    $self->status_message("Got metagenomic_nucleotide_model ".$metagenomic_nucleotide_model->__display_name__) if $metagenomic_nucleotide_model;
    my $metagenomic_protein_model = $model->metagenomic_protein_model;
    $self->status_message("Got metagenomic_protein_model ".$metagenomic_protein_model->__display_name__) if $metagenomic_protein_model;
    my $viral_nucleotide_model = $model->viral_nucleotide_model;
    $self->status_message("Got viral_nucleotide_model ".$viral_nucleotide_model->__display_name__) if $viral_nucleotide_model;
    my $viral_protein_model = $model->viral_protein_model;
    $self->status_message("Got viral_protein_model ".$viral_protein_model->__display_name__) if $viral_protein_model;

    my $screen_contamination = Genome::Model::MetagenomicShotgun::Build::ScreenContamination->create(
        build => $build,
    );
    return if not $screen_contamination;
    return if not $screen_contamination->execute;

    my %mg_nucleotide_results = $self->_run_meta_nt($build, $screen_contamination->unaligned);
    return if not %mg_nucleotide_results;

    my @mg_nucleotide_unaligned = @{$mg_nucleotide_results{unaligned}};
    my @mg_protein_aligned = $self->_run_meta_nr($build,  @mg_nucleotide_unaligned);
    return if not @mg_protein_aligned;

    my @mg_nucleotide_aligned = @{$mg_nucleotide_results{aligned}};
    my $viral_nr_ok = $self->_run_viral_nr($build, @mg_nucleotide_aligned, @mg_protein_aligned);
    return if not $viral_nr_ok;

    my $viral_nt_ok = $self->_run_viral_nt($build, @mg_nucleotide_aligned, @mg_protein_aligned);
    return if not $viral_nt_ok;

    return 1;
}

sub _run_meta_nt {
    my ($self, $build, @cs_unaligned) = @_;

    my $metagenomic_nucleotide_model = $build->model->metagenomic_nucleotide_model;
    my $mg_nucleotide_build = $self->_start_build($metagenomic_nucleotide_model, @cs_unaligned);
    my $mg_nt_build_ok = $self->_wait_for_build($mg_nucleotide_build);
    return if not $mg_nt_build_ok;
    my $link_alignments = $self->_link_sub_build_alignments_to_build(build => $build, sub_build => $mg_nucleotide_build, sub_model_name => 'metagenomic_nucleotide');
    return if not $link_alignments;

    my @mg_nucleotide_aligned = $self->extract_data($mg_nucleotide_build, "aligned");
    my @mg_nucleotide_unaligned = $self->_extract_data($mg_nucleotide_build, "unaligned");

    return ( 
        aligned => \@mg_nucleotide_aligned,
        unaligned => \@mg_nucleotide_unaligned
    );
}

sub _run_meta_nr {
    my ($self, $build, @mg_nucleotide_unaligned) = @_;

    my $metagenomic_protein_model = $build->model->metagenomic_protein_model;
    my $mg_protein_build = $self->_start_build($metagenomic_protein_model, @mg_nucleotide_unaligned);
    my $mg_nr_build_ok = $self->_wait_for_build($mg_protein_build);
    return if not $mg_nr_build_ok;

    my $link_alignments = $self->_link_sub_build_alignments_to_build(build => $build, sub_build => $mg_protein_build, sub_model_name => 'metagenomic_protein');
    return if not $link_alignments;

    my @mg_protein_aligned = $self->_extract_data->($mg_protein_build, "aligned");

    return @mg_protein_aligned;
}

sub _run_viral_nt {
    my ($self, $build, @instrument_data);

    my $viral_nucleotide_model = $build->model->viral_nucleotide_model;
    my $viral_nucleotide_build = $self->_start_build($viral_nucleotide_model, @instrument_data);
    my $viral_nt_build_ok = $self->_wait_for_build($viral_nucleotide_build);
    return if not $viral_nt_build_ok;

    my $link_alignments = $self->_link_sub_build_alignments_to_build(build => $build, sub_build => $viral_nucleotide_build, sub_model_name => 'viral_nucleotide');
    return if not $link_alignments;

    return 1;
}

sub _run_viral_nr {
    my ($self, $build, @instrument_data) = @_;

    my $viral_protein_model = $build->model->viral_protein_model;
    my $viral_protein_build = $self->_start_build($viral_protein_model, @instrument_data);
    my $viral_nr_build_ok = $self->_wait_for_build($viral_protein_build);
    return if not $viral_nr_build_ok;

    my $link_alignments = $self->_link_sub_build_alignments_to_build(build => $build, sub_build => $viral_protein_build, sub_model_name => 'viral_protein');
    return if not $link_alignments;

    return 1;
}

1;

