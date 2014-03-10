package Genome::Model::MetagenomicCompositionShotgun::Command::Status;

use Genome;
use strict;
use warnings;

class Genome::Model::MetagenomicCompositionShotgun::Command::Status {
    is => 'Genome::Model::MetagenomicCompositionShotgun::Command',
    doc => 'display status of sub-models of a metagenomic shotgun composition build',
    has => [
        item => {
            is => 'Text',
            doc => 'name or id of build or model',
            is_optional => 1,
            shell_args_position => 1,
        },
    ],
};

sub execute {
    my $self = shift;
    unless ($self->item) {
        $self->error_message("Please specify an item.");
        exit;
    }

    my ($mcs_model, $mcs_build);
    # If the ITEM is a number then try to get the model or build by ID.
    if ($self->item =~ /^\d+$/) {
        $mcs_model = Genome::Model->get($self->item);
        if ($mcs_model) {
            $self->debug_message("Found model " . $mcs_model->name . " using model ID " . $self->item . ".");
        }
        else {
            $mcs_build = Genome::Model::Build->get($self->item);
            $mcs_model = $mcs_build->model;
            $self->debug_message("Found model " . $mcs_model->name . " (" . $mcs_model->id . ") from build ID " . $self->item . ".");
        }
    }
    # If the ITEM is not a number or it failed to get model/build using a number ID then try as name.
    if (!$mcs_model && ($mcs_model = Genome::Model->get(name => $self->item))) {
        $self->debug_message("Found model using Name=" . $self->item . "; model ID is " . $mcs_model->id . ".");
    }
    # If we got the model from ITEM then get the latest build.
    if (!$mcs_build && $mcs_model) {
        $mcs_build = $mcs_model->latest_build;
        $self->debug_message("Using latest build from model " . $mcs_model->name . "; build ID is " . $mcs_build->id . ".");
    }
    unless ($mcs_build && $mcs_model) {
        $self->error_message("Failed to get build and model.");
        exit;
    }

    my $hcs_model = $mcs_model->_contamination_screen_alignment_model;
    my @meta_models = $mcs_model->_metagenomic_alignment_models;

    my $hcs_build = $hcs_model->latest_build;
    my $hcs_status = $hcs_build->status if ($hcs_build);
    $hcs_status = 'Not running' unless($hcs_status);
    $self->debug_message($hcs_model->name . ": " . $hcs_build->status) if ($hcs_build);
    for my $meta_model (@meta_models) {
        my $meta_build = $meta_model->latest_build;
        my $meta_status = $meta_build->status if ($meta_build);
        $meta_status = 'Not running' unless($meta_status);
        $self->debug_message($meta_model->name . ": " . $meta_status);
    }
    $self->debug_message($mcs_model->name . ": " . $mcs_build->status);

    return 1;
}

1;
