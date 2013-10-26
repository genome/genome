package Genome::Model::SomaticValidation::Command::Loh;

use strict;
use warnings;

use Genome;

class Genome::Model::SomaticValidation::Command::Loh {
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

    return 1 if $self->should_skip_run;

    return;

    #ideally would try to shortcut...
    #...but need to run dispatcher to get control result

    #my @params = $self->_params_for_result;
    #return unless @params;

    #my $result = Genome::Model::Tools::DetectVariants2::Classify::Loh->get_with_lock(@params);
    #return unless $result;

    #$self->status_message('Using existing result: ' . $result->id);

    #return $self->link_result_to_build($result);
}

sub sub_command_category { 'pipeline steps' }

sub execute {
    my $self = shift;
    my $build = $self->build;

    return 1 if $self->should_skip_run;

    my @params = $self->_params_for_result;
    if(@params) {
        my $result = Genome::Model::Tools::DetectVariants2::Classify::Loh->get_or_create(@params);
        return $self->link_result_to_build($result);
    } else {
        die $self->error_message('Failed to gather params to generate result');
    }

    $self->status_message("Identify LOH step completed");
    return 1;
}

sub should_skip_run {
    my $self = shift;
    my $build = $self->build;

    unless ($build){
        die $self->error_message("no build provided!");
    }

    unless(defined($build->model->loh_version)){
        $self->status_message("No LOH version was found, skipping LOH detection!");
        return 1;
    }

    unless(defined($build->model->snv_detection_strategy)){
        $self->status_message("No SNV Detection Strategy, skipping LOH.");
        return 1;
    }

    unless(defined($build->model->loh_snv_detection_strategy)) {
        $self->status_message("No LOH SNV detection strategy, skipping LOH");
        return 1;
    }

    unless($build->normal_sample) {
        $self->status_message('No control to compare against, skipping LOH');
        return 1;
    }

    return;
}

sub _fetch_control_snv_result {
    my $self = shift;
    my $build = $self->build;

    my $existing_result_user = Genome::SoftwareResult::User->get(user => $build, label => 'control_snv_result_for_loh');
    if($existing_result_user) {
        return $existing_result_user->software_result;
    }

    return unless $build->model->loh_snv_detection_strategy;

    my $alignment_result = $build->control_merged_alignment_result;
    my $bam = $alignment_result->merged_alignment_bam_path;

    my $dir = $build->data_directory .'/control_variants_for_loh';
    Genome::Sys->create_directory($dir);

    my %params = (
        snv_detection_strategy => $build->model->loh_snv_detection_strategy,
        aligned_reads_input => $bam,
        aligned_reads_sample => $build->normal_sample->name,
        reference_build_id => $build->reference_sequence_build->id,
        output_directory => $dir,
    );

    my $command = Genome::Model::Tools::DetectVariants2::Dispatcher->create(%params);
    my $rv = $command->execute;
    my $err = $@;
    unless ($rv){
        die $self->error_message("Failed to execute detect variants dispatcher(err:$@) with params:\n".Data::Dumper::Dumper \%params);
    }

    my $result = $command->snv_result;
    unless(grep($_ eq $result, $build->results)) {
        $result->add_user(user => $build, label => 'control_snv_result_for_loh');
    }

    return $result;
}


sub _params_for_result {
    my $self = shift;
    my $build = $self->build;

    my $prior_result_user = Genome::SoftwareResult::User->get(user => $build, label => 'snv_result');
    my $prior_result = $prior_result_user->software_result;
    my $control_result = $self->_fetch_control_snv_result;
    my $loh_version = $build->model->loh_version;

    return unless $prior_result and $control_result;

    return (
        prior_result_id => $prior_result->id,
        control_result_id => $control_result->id,
        classifier_version => $loh_version,
        test_name => $ENV{GENOME_SOFTWARE_RESULT_TEST_NAME} || undef,
    );
}

sub link_result_to_build {
    my $self = shift;
    my $result = shift;
    my $build = $self->build;

    $result->add_user(user => $build, label => 'loh');
    Genome::Sys->create_symlink($result->output_dir, join('/', $build->data_directory, 'loh'));

    return 1;
}

1;

