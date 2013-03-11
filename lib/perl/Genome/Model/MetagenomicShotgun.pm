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

    my @original_instrument_data = $build->instrument_data;
    my $screen_contamination = Genome::Model::MetagenomicShotgun::Build::AlignTo->create(
        sub_model_label => 'contamination_screen_model',
        build => $build,
        instrument_data => \@original_instrument_data,
    );
    return if not $screen_contamination;
    return if not $screen_contamination->execute;

    my $meta_nt = Genome::Model::MetagenomicShotgun::Build::AlignTo->create(
        sub_model_label => 'metagenomic_nucleotide_model',
        build => $build,
        instrument_data => [ $screen_contamination->unaligned ],
    );
    return if not $meta_nt;
    return if not $meta_nt->execute;

    my $meta_nr = Genome::Model::MetagenomicShotgun::Build::AlignTo->create(
        sub_model_label => 'metagenomic_protein_model',
        build => $build,
        instrument_data => [ $meta_nt->unaligned ],
    );
    return if not $meta_nr;
    return if not $meta_nr->execute;

    my $viral_nr = Genome::Model::MetagenomicShotgun::Build::AlignTo->create(
        sub_model_label => 'viral_protein_model',
        build => $build,
        instrument_data => [ $meta_nt->aligned, $meta_nr->aligned ],
    );
    return if not $viral_nr;
    return if not $viral_nr->execute;

    my $viral_nt = Genome::Model::MetagenomicShotgun::Build::AlignTo->create(
        sub_model_label => 'viral_nucleotide_model',
        build => $build,
        instrument_data => [ $meta_nt->aligned, $meta_nr->aligned ],
    );
    return if not $viral_nt;
    return if not $viral_nt->execute;

    return 1;
}

1;

