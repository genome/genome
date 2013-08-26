package Genome::Model::SomaticValidation::Command::IdentifyPreviouslyDiscoveredVariations;

use strict;
use warnings;

use Genome;

class Genome::Model::SomaticValidation::Command::IdentifyPreviouslyDiscoveredVariations{
    is => 'Genome::Command::Base',
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
        }
    ],
    has_param => [
        lsf_queue => {
            default => 'apipe',
        },
    ],
};

sub shortcut {
    my $self = shift;

    if($self->should_skip_run) {
        return 1;
    }

    my $build = $self->build;

    my @results;
    my $expected_result_count = 0;

    for my $type ('snv', 'indel') {
        my $accessor = join('_', $type, 'detection_strategy');
        if($build->model->$accessor) {
            $expected_result_count++;
            my @params_for_result = $self->params_for_result($type);
            unless(@params_for_result) {
                die $self->error_message('Could not get params for ' . $type);
            }
            my $result = Genome::Model::Tools::DetectVariants2::Classify::PreviouslyDiscovered->get_with_lock(@params_for_result);
            push @results, $result if $result;
        }
    }

    unless($expected_result_count > 0 && scalar(@results) == $expected_result_count) {
        return;
    }

    for my $r (@results) {
        $self->link_result_to_build(@results) or die('Failed to link result.');
    }

    return 1;
}

sub sub_command_category { 'pipeline steps' }

sub execute {
    my $self  = shift;

    if($self->should_skip_run) {
        return 1;
    }

    my $build = $self->build;

    $self->status_message("Comparing detected variants to previously discovered variations");

    my ($snv_result, $indel_result);

    my $prev_variations_build = $build->previously_discovered_variations_build;
    $snv_result   = $prev_variations_build->snv_result;
    $indel_result = $prev_variations_build->indel_result;

    unless ($indel_result or $snv_result) {
        die $self->error_message("No indel or snv result found on previously discovered variations build. This is unsupported!  Failing.");
    }

    my $version = 2;
    #my $version = GMT:BED:CONVERT::version();  TODO, something like this instead of hardcoding

    if ($build->snv_detection_strategy){
        if($snv_result) {
            my @params = $self->params_for_result('snv');
            my $result = Genome::Model::Tools::DetectVariants2::Classify::PreviouslyDiscovered->get_or_create(@params);
            unless($result) {
                die $self->error_message('Failed to generate result for SNVs.');
            }
            $self->link_result_to_build($result);
        } else {
            $self->status_message("No snv result found on previously discovered variations build, skipping snv intersection");
        }
    }

    if ($build->indel_detection_strategy) {
        if($indel_result) {
            my @params = $self->params_for_result('indel');
            my $result = Genome::Model::Tools::DetectVariants2::Classify::PreviouslyDiscovered->get_or_create(@params);
            unless($result) {
                die $self->error_message('Failed to generate result for indels.');
            }
            $self->link_result_to_build($result);
        } else {
            $self->status_message("No indel result found on previously discovered variations build, skipping indel intersection");
        }
    }

    $self->status_message("Identify Previously Discovered Variations step completed");
    return 1;
}

sub should_skip_run {
    my $self = shift;
    my $build = $self->build;

    unless ($build){
        die $self->error_message("no build provided!");
    }

    unless(defined($build->model->snv_detection_strategy) or defined($build->model->indel_detection_strategy)){
        $self->status_message("No SNV or indel detection strategy, skipping identify previously discovered variants.");
        return 1;
    }

    my $prev_variations_build = $build->previously_discovered_variations_build;
    unless ($prev_variations_build) {
        $self->warning_message('no previously_discovered_variations_build provided !');
        return 1;
    }

    return;
}

sub params_for_result {
    my $self = shift;
    my $variant_type = shift;
    my $build = $self->build;

    my $label;
    if($variant_type eq 'snv' and $build->loh_version and $build->normal_sample) {
        $label = 'loh';
    } else {
        $label = $variant_type . '_result';
    }
    my $prior_result_user = Genome::SoftwareResult::User->get(label => $label, user => $build);
    my $prior_result = $prior_result_user->software_result;

    my $variations_build = $build->previously_discovered_variations_build;
    my $accessor = $variant_type . '_result';
    my $previously_discovered_result = $variations_build->$accessor;

    unless($prior_result and $previously_discovered_result) {
        return; #can't create a result for old things
    }

    return (
        prior_result_id => $prior_result->id,
        previously_discovered_result_id => $previously_discovered_result->id,
        classifier_version => 1,
        variant_type => $variant_type,
        test_name => $ENV{GENOME_SOFTWARE_RESULT_TEST_NAME} || undef,
        skip_filtering => undef,
    );
}

sub link_result_to_build {
    my $self = shift;
    my $result = shift;
    my $build = $self->build;

    $self->status_message('Linking result ' . $result->id . ' to build.');
    $result->add_user(user => $build, label => $result->variant_type . '_identify_previously_discovered_result');
    Genome::Sys->create_directory($build->data_directory."/novel");

    my $type = $result->variant_type;
    my $novel = $type . 's.hq.novel.v2.bed';
    my $previously_detected = $type . 's.hq.previously_detected.v2.bed';
    my $novel_detected_path      = join('/', $build->data_directory, 'novel', $novel);
    my $previously_detected_path = join('/', $build->data_directory, 'novel', $previously_detected);
    Genome::Sys->create_symlink($result->path($novel), $novel_detected_path);
    Genome::Sys->create_symlink($result->path($previously_detected), $previously_detected_path);

    return 1;
}

1;

