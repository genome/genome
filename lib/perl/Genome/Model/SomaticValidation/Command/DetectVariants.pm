package Genome::Model::SomaticValidation::Command::DetectVariants;

use strict;
use warnings;

use Genome;

class Genome::Model::SomaticValidation::Command::DetectVariants{
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
        },
        bam_path => {
            is => 'Text',
            is_input => 1,
            doc => 'path to tumor bam file',
        },
        control_bam_path => {
            is => 'Text',
            is_input => 1,
            doc => 'path to normal bam file',
            is_optional => 1,
        },
    ],
    has_param => [
        lsf_queue => {
            default => $ENV{GENOME_LSF_QUEUE_BUILD_WORKER},
        },
    ],
    #specific things used by the final process-validation script
    has_output_optional => [
        hq_snv_file => {
            is => 'Text',
            doc => 'The filtered snv file from the DV run',
        },
        lq_snv_file => {
            is => 'Text',
            doc => 'The unfiltered snv file from the DV run',
        }
    ],
    doc => 'detect variants',
};

sub sub_command_category { 'pipeline steps' }

sub execute{
    my $self = shift;

    $self->debug_message("Executing detect variants step");
    my $build = $self->build;
    unless ($build){
        die $self->error_message("no build provided!");
    }

    unless ($build->snv_detection_strategy or $build->indel_detection_strategy or $build->sv_detection_strategy or $build->cnv_detection_strategy) {
        $self->warning_message("No detection strategies provided. Skipping detect variants step");
        return 1;
    }

    my %params;
    $params{snv_detection_strategy} = $build->snv_detection_strategy if $build->snv_detection_strategy;
    $params{indel_detection_strategy} = $build->indel_detection_strategy if $build->indel_detection_strategy;
    $params{sv_detection_strategy} = $build->sv_detection_strategy if $build->sv_detection_strategy;
    $params{cnv_detection_strategy} = $build->cnv_detection_strategy if $build->cnv_detection_strategy;

    my $tumor_bam = $self->bam_path;
    unless (-e $tumor_bam){
        die $self->error_message("No tumor bam found for somatic model");
    }
    $params{aligned_reads_input} = $tumor_bam;

    my $reference_build = $build->reference_sequence_build;
    my $reference_fasta = $reference_build->full_consensus_path('fa');
    unless(-e $reference_fasta){
        die $self->error_message("fasta file for reference build doesn't exist!");
    }
    $params{reference_build_id} = $reference_build->id;

    my $output_dir = $build->data_directory."/variants";
    $params{output_directory} = $output_dir;

    my $aligned_reads_sample = $build->tumor_sample->name;
    $params{aligned_reads_sample} = $aligned_reads_sample;

    if($build->normal_sample) {
        my $normal_bam = $self->control_bam_path;
        unless (-e $normal_bam){
            die $self->error_message("No normal bam found for somatic model");
        }
        $params{control_aligned_reads_input} = $normal_bam;

        my $control_aligned_reads_sample = $build->normal_sample->name;
        $params{control_aligned_reads_sample} = $control_aligned_reads_sample;
    }

    my $command = Genome::Model::Tools::DetectVariants2::Dispatcher->create(%params);
    unless ($command){
        die $self->error_message("Couldn't create detect variants dispatcher from params:\n".Data::Dumper::Dumper \%params);
    }
    my $rv = $command->execute;
    my $err = $@;
    unless ($rv){
        die $self->error_message("Failed to execute detect variants dispatcher(err:$@) with params:\n".Data::Dumper::Dumper \%params);
    }
    else {
        #users set below
    }

    $self->debug_message("detect variants command completed successfully");

    my $version = 2;
    #my $version = GMT:BED:CONVERT::version();  TODO, something like this instead of hardcoding

    if ($build->snv_detection_strategy){
        my $snv_result = $command->snv_result;
        $snv_result->add_user(user => $build, label => 'snv_result');
        my $result = $build->data_set_path("variants/snvs.hq",$version,'bed');
        unless (-e $result){
            my $unexpected_format_output = $command->_snv_hq_output_file;
            unless (-e $unexpected_format_output){
                die $self->error_message("Expected hq detected snvs file $result, but it does not exist!");
            }
            symlink($unexpected_format_output, $result);
        }
        my $lq_snv_result = $command->snv_lq_result;
        if($lq_snv_result) {
            $lq_snv_result->add_user(user => $build, label => 'snv_lq_result');
        }
        my $lq_result = $build->data_set_path("variants/snvs.lq",$version,'bed');
        unless (-e $lq_result){
            #these are not set on the detect variants command, so I'm hardcoding them.
            my $unexpected_filename_output = $build->data_directory."/variants/snvs.lq.bed";
            unless (-e $unexpected_filename_output){
                die $self->error_message("Expected lq detected snvs $unexpected_filename_output, but it does not exist");
            }
            symlink($unexpected_filename_output, $lq_result);
        }
    }

    if ($build->indel_detection_strategy){
        my $indel_result = $command->indel_result;
        $indel_result->add_user(user => $build, label => 'indel_result');
        my $result = $build->data_set_path("variants/indels.hq",$version,'bed');
        unless (-e $result){
            my $unexpected_filename_output = $command->_indel_hq_output_file;
            unless (-e $unexpected_filename_output){
                die $self->error_message("Expected hq detected indels file $result, but it does not exist!");
            }
            symlink($unexpected_filename_output, $result);
        }
        my $lq_indel_result = $command->indel_lq_result;
        if($lq_indel_result) {
            $lq_indel_result->add_user(user => $build, label => 'indel_lq_result');
        }
        my $lq_result = $build->data_set_path("variants/indels.lq",$version,'bed');
        unless (-e $lq_result){
            #these are not set on the detect variants command, so I'm hardcoding them.
            my $unexpected_filename_output = $build->data_directory."/variants/indels.lq.bed";
            if (-e $unexpected_filename_output){
                symlink($unexpected_filename_output, $lq_result);
            } else {
                $self->debug_message("No lq indel file found. Creating an empty file at $lq_result");
                system("touch $lq_result");
            }
        }
    }

    if ($build->sv_detection_strategy){
        my $sv_result = $command->sv_result;
        $sv_result->add_user(user => $build, label => 'sv_result');
        my $result = $build->data_set_path("variants/svs.hq",$version,'bed');
        unless (-e $result){
            my $unexpected_filename_output = $command->_sv_hq_output_file;
            unless (-e $unexpected_filename_output){
                die $self->error_message("Expected hq detected snvs file $result, but it does not exist!");
            }
            symlink($unexpected_filename_output, $result);
        }
    }

    if ($build->cnv_detection_strategy){
        my $cnv_result = $command->cnv_result;
        $cnv_result->add_user(user => $build, label => 'cnv_result');
        my $result = $build->data_set_path("variants/cnvs.hq",$version,'bed');
        unless (-e $result){
            my $unexpected_filename_output = $command->_cnv_hq_output_file;
            unless (-e $unexpected_filename_output){
                die $self->error_message("Expected hq detected snvs file $result, but it does not exist!");
            }
            symlink($unexpected_filename_output, $result);
        }
    }

    #find the relevant output files for outputs from this step
    my $hq_snv_file = $self->build->data_set_path('variants/snvs.hq', $version, 'bed');
    my $lq_snv_file = $self->build->data_set_path('variants/snvs.lq', $version, 'bed');

    unless($hq_snv_file) {
        die $self->error_message('Could not find an HQ snv file.');
    }
    unless($lq_snv_file) {
        die $self->error_message('Could not find an LQ snv file.');
    }

    $self->hq_snv_file($hq_snv_file);
    $self->lq_snv_file($lq_snv_file);

    $self->debug_message("detect variants step completed");

    return 1;
}

1;
