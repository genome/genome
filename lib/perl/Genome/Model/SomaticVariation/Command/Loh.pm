package Genome::Model::SomaticVariation::Command::Loh;

use strict;
use warnings;
use Genome;
use Genome::Info::IUB;
use Genome::Model::Tools::DetectVariants2::Utilities qw(
    final_result_for_variant_type
);

class Genome::Model::SomaticVariation::Command::Loh {
    is => 'Genome::Command::Base',
    has =>[
        build_id => {
            is => 'Text',
            is_input => 1,
            is_output => 1,
            doc => 'build id of SomaticVariation model',
        },
        build => {
            is => 'Genome::Model::Build::SomaticVariation',
            id_by => 'build_id',
        }
    ],
    has_optional => [
        variant_bed_file => {
            is => 'Text',
            is_input => 1,
            doc => 'use this to specify a specific file as an input',
        },
        output_directory => {
            is => 'Text',
            is_input => 1,
            doc => 'use this to send the output of LOH to an alternate directory',
        },
    ],
    has_param => [
        lsf_queue => {
            default => $ENV{GENOME_LSF_QUEUE_BUILD_WORKER_ALT},
        },
    ],
};

sub shortcut {
    my $self = shift;

    return 1 if $self->should_skip_run;

    my @params = $self->_params_for_result;
    return unless @params;

    my $result = Genome::Model::Tools::DetectVariants2::Classify::Loh->get_with_lock(@params);
    return unless $result;

    $self->debug_message('Using existing result: ' . $result->id);

    return $self->link_result_to_build($result);
}

sub execute {
    my $self = shift;
    my $build = $self->build;

    return 1 if $self->should_skip_run;

    my $normal_build = $build->normal_build;
    unless ($normal_build){
        die $self->error_message("No previous normal build found on somatic build!");
    }

    my @params = $self->_params_for_result;
    if(@params) {
        my $result = Genome::Model::Tools::DetectVariants2::Classify::Loh->get_or_create(@params);
        return $self->link_result_to_build($result);
    }
    #else using old data...fall back on running by hand

    #Use snvs bed for intersecting with bed style outputs from detect-variants
    my $normal_snvs;
    # Wrap this call in an eval so the whole process won't die if there's no bed, 
    # we will fall back on the annotation format, and temporarily convert it
    eval {
        $normal_snvs = $normal_build->filtered_snvs_bed; 
    };
    unless(defined($normal_snvs) && -e $normal_snvs){
        # in the event that no filtered_snvs_bed is found, fall back on samtools output, and convert that to bed
        unless($normal_snvs = $self->get_temp_bed_snvs($normal_build)){
            die $self->error_message("Could not find snvs_bed from normal reference-alignment build.");
        }
    }
    $self->debug_message("Looking for LOH events in SNV output");

    my $version = 2;

    my $detected_snv_path = defined($self->variant_bed_file) ? $self->variant_bed_file : $build->data_set_path("variants/snvs.hq",$version,'bed'); 
    my $output_dir = defined($self->output_directory) ? $self->output_directory : $build->data_directory."/loh" ;
    unless(Genome::Sys->create_directory($output_dir)){
        die $self->error_message("Failed to create the ./loh subdir");
    }

    my $somatic_output = $output_dir."/snvs.somatic.v".$version.".bed";
    my $loh_output = $output_dir."/snvs.loh.v".$version.".bed";

    my $aligned_reads_input = $build->tumor_build->whole_rmdup_bam_file;
    my $control_aligned_reads_input = $build->normal_build->whole_rmdup_bam_file;
    my $reference_build_id = $build->reference_sequence_build->id;

    Genome::Model::Tools::DetectVariants2::Classify::Loh->run_loh( $normal_snvs, $detected_snv_path, $somatic_output, $loh_output );

    $self->debug_message("Identify LOH step completed");
    return 1;
}

sub should_skip_run {
    my $self = shift;
    my $build = $self->build;

    unless ($build){
        die $self->error_message("no build provided!");
    }

    unless(defined($build->loh_version)){
        $self->debug_message("No LOH version was found, skipping LOH detection!");
        return 1;
    }

    unless(defined($build->model->snv_detection_strategy)){
        $self->debug_message("No SNV Detection Strategy, skipping LOH.");
        return 1;
    }

    return;
}

# return a path to a temp file containing a bed version of the samtools style snv_file
sub get_temp_bed_snvs {
    my $self = shift;
    my $normal_build = shift;
    my $normal_snvs = $normal_build->snv_file;
    unless (-e $normal_snvs) {
        return;
    }

    my $temp_bed_file = Genome::Sys->create_temp_file_path;
    my $convert = Genome::Model::Tools::Bed::Convert::Snv::SamtoolsToBed->create( 
                        source => $normal_snvs, 
                        output => $temp_bed_file);

    unless($convert->execute){
        die $self->error_message("Failed to run conversion from samtools to bed format");
    }
    return $temp_bed_file;
}

sub _params_for_result {
    my $self = shift;
    my $build = $self->build;

    my $prior_result = final_result_for_variant_type([$build->results], 'snv');
    my $control_result = final_result_for_variant_type([$build->normal_build->results], 'snv');
    my $loh_version = $build->loh_version;

    return unless $prior_result and $control_result; #can't create a result for old things

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

    # If an output directory has been provided (not going to the build dir), the build is not a user
    my $symlink_target;
    if ($self->output_directory) {
        $symlink_target = $self->output_directory;
    } else {
        my $build = $self->build;
        $result->add_user(user => $build, label => 'uses');
        $symlink_target = join('/', $build->data_directory, 'loh')
    }

    Genome::Sys->create_symlink($result->output_dir, $symlink_target);

    return 1;
}

1;

