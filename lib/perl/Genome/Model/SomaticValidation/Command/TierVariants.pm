package Genome::Model::SomaticValidation::Command::TierVariants;

use strict;
use warnings;

use Genome;

class Genome::Model::SomaticValidation::Command::TierVariants{
    is => 'Command::V2',
    has =>[
        build_id => {
            is => 'Text',
            is_input => 1,
            is_output => 1,
            doc => 'build id of SomaticValidation model',
        },
        build => {
            is => 'Genome::Model::Build::SomaticValidation',
            id_by => 'build_id',
        },
    ],
    has_param => [
        lsf_queue => {
            default => 'apipe',
        },
    ],
};

sub shortcut {
    my $self = shift;
    my $build = $self->build;

    return 1 if $self->should_skip_run;

    my @existing_results;
    for my $qual ('hq', 'lq') {
        TYPE: for my $variant_type ('snv', 'indel') {
            my $accessor = $variant_type . '_detection_strategy';
            next TYPE unless $build->$accessor;

            my @params = $self->params_for_result($qual, $variant_type);
            unless(@params) {
                die $self->error_message('Failed to get params for ' . $qual . ' ' . $variant_type . ' tiering');
            }

            my $existing_result = Genome::Model::Tools::DetectVariants2::Classify::Tier->get_with_lock(@params);
            if($existing_result) {
                push @existing_results, $existing_result;
            } else {
                return;
            }
        }
    }

    for my $existing_result (@existing_results) {
        $self->link_result_to_build($existing_result);
    }

    return 1;
}

sub should_skip_run {
    my $self = shift;
    my $build = $self->build;

    unless($self->build) {
        die $self->error_message('no build provided');
    }

    unless($build->tiering_version) {
        $self->status_message('No tiering version specified... skipping run.');
        return 1;
    }

    unless($build->annotation_build) {
        $self->status_message('No annotation build specified... skipping run.');
        return 1;
    }

    return;
}

sub sub_command_category { 'pipeline steps' }

sub execute {
    my $self = shift;
    my $build = $self->build;
    unless ($build){
        die $self->error_message("no build provided!");
    }

    return 1 if $self->should_skip_run;

    $self->status_message("executing tier variants step on snvs and indels");

    for my $qual ('hq', 'lq') {
        TYPE: for my $variant_type ('snv', 'indel') {
            my $accessor = $variant_type . '_detection_strategy';
            next TYPE unless $build->$accessor;

            my @params = $self->params_for_result($qual, $variant_type);
            if(@params) {
                my $result = Genome::Model::Tools::DetectVariants2::Classify::Tier->get_or_create(@params);
                if($result) {
                    $self->link_result_to_build($result);
                } else {
                    die $self->error_message('Failed to generate result for ' . $qual . ' ' . $variant_type);
                }
            } else {
                unless(@params) {
                    die $self->error_message('Failed to get params for ' . $qual . ' ' . $variant_type . ' tiering');
                }
            }
        }
    }

    $self->status_message("Tier Variants step completed");
    return 1;
}

sub params_for_result {
    my $self = shift;
    my $qual = shift;
    my $variant_type = shift;
    my $build = $self->build;

    my $label;
    if($qual eq 'hq') {
        my $accessor = $variant_type . '_result';
        if($build->previously_discovered_variations_build and $build->previously_discovered_variations_build->$accessor) {
            $label = $variant_type . '_identify_previously_discovered_result';
        } elsif($variant_type eq 'snv' and $build->model->loh_version) {
            $label = 'loh';
        } else {
            $label = $variant_type . '_result';
        }
    } elsif($qual eq 'lq') {
        $label = $variant_type . '_lq_result';
    } else {
        die('Unknown quality passed to params_for_result: ' . $qual);
    }

    my $user = Genome::SoftwareResult::User->get(user => $build, label => $label);
    unless($user) {
        $self->error_message('No result matching label ' . $label);
        return;
    }
    my $result = $user->software_result;
    return unless $result;

    return (
        variant_type => $variant_type,
        prior_result_id => $result->id,
        annotation_build_id => $build->annotation_build->id,
        classifier_version => $build->tiering_version,
        test_name => $ENV{GENOME_SOFTWARE_RESULT_TEST_NAME} || undef,
    );
}

sub link_result_to_build {
    my $self = shift;
    my $result = shift;
    my $build = $self->build;

    my $type = 'hq';
    my $opposite_type = 'lq';
    my $prior_result = $result->prior_result;
    if($prior_result->isa('Genome::Model::Tools::DetectVariants2::Result::Combine::LqUnion')) {
        $type = 'lq';
        $opposite_type = 'hq';
    }

    $result->add_user(label => $result->variant_type . '_tiered_' . $type, user => $build);

    my $effects_dir = $build->data_directory."/effects";
    unless(-d $effects_dir){
        Genome::Sys->create_directory($effects_dir);
    }

    for my $f (glob($result->output_dir . '/*')) {
        my $name = File::Basename::fileparse($f);
        next unless $name !~ $opposite_type;

        Genome::Sys->create_symlink($f, join('/', $effects_dir, $name));
    }

    return 1;
}

1;

