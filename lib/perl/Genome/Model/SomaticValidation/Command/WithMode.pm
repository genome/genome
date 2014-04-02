package Genome::Model::SomaticValidation::Command::WithMode;

use strict;
use warnings;
use Genome;

class Genome::Model::SomaticValidation::Command::WithMode {
    is => 'Command::V2',
    has_input => [
        build_id => {
            is => 'Number',
            doc => 'ID of the SomaticValidation build upon which to run coverage stats',
        },
        mode => {
            is => 'Text',
            valid_values => ['tumor', 'normal'],
        },
    ],
    has => [
        build => {
            is => 'Genome::Model::Build::SomaticValidation',
            id_by => 'build_id',
            doc => 'The build upon which to run coverage stats',
        },
    ],
};

sub should_run {
    my $self = shift;

    return $self->sample_for_mode;
}

sub sample_for_mode {
    my $self = shift;
    my $mode = $self->mode;
    my $sample_acc = $mode . '_sample';
    return $self->build->$sample_acc;
}

sub alignment_result_for_mode {
    my $self = shift;
    my $alignment_result;
    if($self->mode eq 'tumor') {
        $alignment_result = $self->build->merged_alignment_result;
    } else {
        $alignment_result = $self->build->control_merged_alignment_result;
    }

    unless($alignment_result) {
        die $self->error_message('No alignment result found for ' . $self->mode);
    }
    return $alignment_result;
}

sub link_result_to_build {
    my $self = shift;
    my $result = shift;
    my $subdir_name = shift;
    my $label_base = shift;
    my $build = $self->build;

    my $base_dir = File::Spec->join($build->data_directory, $subdir_name);
    unless (-e $base_dir) {
        Genome::Sys->create_directory($base_dir);
    }

    my $link = File::Spec->join($base_dir, $self->mode);
    my $label = join('_', $label_base, $self->mode);
    Genome::Sys->create_symlink($result->output_dir, $link);
    $result->add_user(label => $label, user => $build);
    $self->add_metrics_to_build($result);
    return 1;
}

sub add_metrics_to_build {
    my $self = shift;
    my $result = shift;

    for my $metric (map {$_->metric_name} $result->metrics) {
        $self->build->set_metric($metric."-".$self->mode, $result->$metric);
    }
    return 1;
}

1;

