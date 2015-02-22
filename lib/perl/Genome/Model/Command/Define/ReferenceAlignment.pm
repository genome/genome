package Genome::Model::Command::Define::ReferenceAlignment;

use strict;
use warnings;

use Genome;

require Carp;
use Regexp::Common;

class Genome::Model::Command::Define::ReferenceAlignment {
    is => 'Genome::Model::Command::Define::HelperDeprecated',
    has => [
        reference_sequence_build => {
            is => 'Genome::Model::Build::ImportedReferenceSequence',
            doc => 'ID or name of the reference sequence to align against',
            example_values => [101947881, 'GRCh37-lite-build37'],
            is_input => 1,
        },
    ],
    has_optional => [
        annotation_reference_build => {
            is => 'Genome::Model::Build::ImportedAnnotation',
            doc => 'ID or name of the the build containing the reference transcript set used for variant annotation',
            is_input => 1,
        },
        genotype_microarray_model => {
            is => 'Genome::Model::GenotypeMicroarray',
            doc => 'ID or name of the genotype microarray model which will be used to obtain the gold snp and genotype files',
            is_optional => 1,
        },
        dbsnp_model => {
            is => 'Genome::Model::ImportedVariationList',
            doc => 'ID or name of the dbSNP model to compare against (the latest build will be selected)',
            is_input => 1,
        },
        dbsnp_build => {
            is => 'Genome::Model::Build::ImportedVariationList',
            doc => 'ID or name of the dbSNP build to compare against',
            is_input => 1,
        },
        target_region_set_names => {
            is => 'Text',
            is_many => 1,
            doc => 'limit the model to take specific capture or PCR instrument data',
        },
        region_of_interest_set_name => {
            is => 'Text',
            doc => 'limit coverage and variant detection to within these regions of interest',
        },
        merge_roi_set => {
            is => 'Boolean',
            default_value => 1,
            doc => 'A flag to merge the region_of_interest_set_name BED file before evaluating coverage.',
        },
        short_roi_names => {
            is => 'Boolean',
            default_value => 1,
            doc => 'A flag to use short ROI names in the BED file used to evaluate coverage.',
        },
        roi_track_name => {
            is => 'Text',
            valid_values => ['target_region','tiled_region'],
            doc => 'If the ROI used is multi-tracked, select one of the tracks for use as the ROI.',
        },
    ],
};

sub resolve_dbsnp {
    my ($self, $rsb) = @_;
    return $self->_resolve_param('dbsnp_build') if $self->dbsnp_build;

    my $dbsnp_model = $self->dbsnp_model;
    if (!$dbsnp_model) {
        $dbsnp_model = Genome::Model::ImportedVariationList->dbsnp_model_for_reference($rsb);
        if (!$dbsnp_model) {
            $self->status_message("no dbsnp_model found.");
            return;
        }
    } else {
        $dbsnp_model = $self->_resolve_param('dbsnp_model');
        if (!defined $dbsnp_model) {
            die $self->error_message("Failed to resolve dbsnp_model identified by " . $self->dbsnp_model);
        }
    }
    my $b = $dbsnp_model->last_complete_build;
    if (!defined $b) {
        die $self->error_message("Failed to find a complete build for dbsnp model " . $dbsnp_model->__display_name__);
    }

    return $b;
}

sub type_specific_parameters_for_create {
    my $self = shift;
    my $rsb = $self->_resolve_param('reference_sequence_build');
    my $arb = $self->_resolve_param('annotation_reference_build');
    my $gmm = $self->_resolve_param('genotype_microarray_model');
    my $dbsnp = $self->resolve_dbsnp($rsb);
    if ($dbsnp && !$rsb->is_compatible_with($dbsnp->model->reference)) {
        die $self->error_message("dbSNP build " . $dbsnp->__display_name__ . " has reference " . $dbsnp->reference->__display_name__ .
            " which does not match the specified reference " . $rsb->__display_name__);
    }
    
    my %params = $self->SUPER::type_specific_parameters_for_create();
    delete $params{dbsnp_model};
    my @params = %params;
    push(@params, reference_sequence_build => $rsb) if $rsb;
    push(@params, annotation_reference_build => $arb) if $arb;
    push(@params, dbsnp_build => $dbsnp) if $dbsnp;
    push(@params, genotype_microarray_model => $gmm) if $gmm;
    return @params;
}

sub listed_params {
    return qw/ id name subject_name subject_type processing_profile_id processing_profile_name reference_sequence_name annotation_reference_name /;
}

sub execute {
    my $self = shift;

    if ($self->dbsnp_build and $self->dbsnp_model and $self->dbsnp_build->model->id != $self->dbsnp_model->id) {
        $self->error_message("Specify one of --dbsnp-build or --dbsnp-model, not both");
        return;
    }
    
    my $result = $self->SUPER::_execute_body(@_);
    return unless $result;

    my $model = Genome::Model->get($self->result_model_id);
    unless ($model) {
        $self->error_message("No model generated for " . $self->result_model_id);
        return;
    }

    # LIMS is preparing actual tables for these in the dw, until then we just manage the names.
    my @target_region_set_names = $self->target_region_set_names;
    if (@target_region_set_names) {
        for my $name (@target_region_set_names) {
            my $i = $model->add_input(value_class_name => 'UR::Value', value_id => $name, name => 'target_region_set_name');
            if ($i) {
                $self->status_message("Modeling instrument-data from target region '$name'");
            }
            else {
                $self->error_message("Failed to add target '$name'!");
                $model->delete;
                return;
            }
        }
    }
    else {
        $self->status_message("Modeling whole-genome (non-targeted) sequence.");
    }
    if ($self->region_of_interest_set_name) {
        my $name = $self->region_of_interest_set_name;
        my $feature_list = Genome::FeatureList->get(name => $name);
        unless ($feature_list) {
            $self->error_message('Failed to find region_of_interest_set_name : '. $name);
            $model->delete;
            return;
        }
        my $i = $model->add_input(value_class_name => 'UR::Value', value_id => $name, name => 'region_of_interest_set_name');
        if ($i) {
            $self->status_message("Analysis limited to region of interest set '$name'");
        }
        else {
            $self->error_message("Failed to add region of interest set '$name'!");
            $model->delete;
            return;
        }
        if (defined($self->merge_roi_set)) {
            my $merge_input = $model->add_input(value_class_name => 'UR::Value', value_id => $self->merge_roi_set, name => 'merge_roi_set');
            if ($merge_input) {
                if ($merge_input->value_id) {
                    $self->status_message('Region of interest set '. $name .' will be merged before analysis.');
                } else {
                    $self->status_message('Region of interest set '. $name .' will NOT be merged before analysis.');
                }
            } else {
                $self->error_message('Failed to set model input merge_roi_set to '. $self->merge_roi_set .'!');
                $model->delete;
                return;
            }
        }
        if (defined($self->short_roi_names)) {
            my $short_names_input = $model->add_input(value_class_name => 'UR::Value', value_id => $self->short_roi_names, name => 'short_roi_names');
            if ($short_names_input) {
                if ($short_names_input->value_id) {
                    $self->status_message('Short names will be used for region of interest set '. $name .'.');
                } else {
                    $self->status_message('The original ROI names will be used for region of interest set '. $name .'.');
                }
            } else {
                $self->error_message('Failed to set model input short_roi_names to '. $self->short_roi_names .'!');
                $model->delete;
                return;
            }
        }
        if (defined($self->roi_track_name)) {
            unless ($feature_list->is_multitracked) {
                $self->error_message('No need to set roi_track_name since region_of_interest_set_name '. $name .' is not multi-tracked!');
                $model->delete;
                return;
            }
            my $roi_track_name_input = $model->add_input(value_class_name => 'UR::Value', value_id => $self->roi_track_name, name => 'roi_track_name');
            if ($roi_track_name_input) {
                if ($roi_track_name_input->value_id) {
                    $self->status_message('Track name '. $roi_track_name_input->value_id .' will be used for region of interest set '. $name .'.');
                }
            } else {
                $self->error_message('Failed to set model input roi_track_name to '. $self->roi_track_name .'!');
                $model->delete;
                return;
            }
        }
    } else {
        $self->status_message("Analyzing whole-genome (non-targeted) reference.");
    }

    return $result;
}

sub _resolve_param {
    my ($self, $param) = @_;

    my $param_meta = $self->__meta__->property($param);
    Carp::confess("Request to resolve unknown property '$param'.") if (!$param_meta);
    my $param_class = $param_meta->data_type;

    my $value = $self->$param;
    return unless $value; # not specified
    if (ref($value) eq 'HASH') {
        # the declaration works with the new UR
        # but this will allow the class to work with the old one also
        # since it did things which are now automatic
        $value = $value->{name};
    }
    return $value if ref($value); # already an object

    my @objs = $self->resolve_param_value_from_text($value, $param_class);
    if (@objs != 1) {
        Carp::confess("Unable to find unique $param_class identified by '$value'. Results were:\n" .
            join('\n', map { $_->__display_name__ . '"' } @objs ));
    }
    $self->$param($objs[0]);
    return $self->$param;
}

1;

