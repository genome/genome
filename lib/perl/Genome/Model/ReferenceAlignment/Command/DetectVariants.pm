package Genome::Model::ReferenceAlignment::Command::DetectVariants;

use strict;
use warnings;

use File::Spec;

use Genome;

class Genome::Model::ReferenceAlignment::Command::DetectVariants {
    is => 'Genome::Model::ReferenceAlignment::Command::PipelineBase',
    has_param => [
        lsf_queue => {
            default => Genome::Config::get('lsf_queue_build_worker'),
        },
        lsf_resource => {
            default => Genome::Config::get('lsf_resource_dv2_dispatcher'),
        },
    ],
};

sub shortcut {
    my $self = shift;

    return $self->_should_skip_run;
}

sub execute {
    my $self = shift;

    return 1 if $self->_should_skip_run;

    $self->debug_message("Executing detect variants step");
    my $build = $self->build;

    my $processing_profile = $build->processing_profile;

    my %params;
    $params{snv_detection_strategy} = $processing_profile->snv_detection_strategy if $processing_profile->snv_detection_strategy;
    $params{indel_detection_strategy} = $processing_profile->indel_detection_strategy if $processing_profile->indel_detection_strategy;
    $params{sv_detection_strategy} = $processing_profile->sv_detection_strategy if $processing_profile->sv_detection_strategy;
    $params{cnv_detection_strategy} = $processing_profile->cnv_detection_strategy if $processing_profile->cnv_detection_strategy;

    my $bam = $build->whole_rmdup_bam_file;
    $params{aligned_reads_input} = $bam;

    my $reference_build = $build->reference_sequence_build;
    my $reference_fasta = $reference_build->full_consensus_path('fa');
    unless(-e $reference_fasta){
        die $self->error_message("fasta file for reference build doesn't exist!");
    }
    $params{reference_build_id} = $reference_build->id;

    my $output_dir = $build->variants_directory;
    $params{output_directory} = $output_dir;

    my $aligned_reads_sample = $build->subject->name;
    $params{aligned_reads_sample} = $aligned_reads_sample;

    $params{result_users} = Genome::SoftwareResult::User->user_hash_for_build($build);

    my $command = Genome::Model::Tools::DetectVariants2::Dispatcher->create(%params);
    unless ($command){
        $self->fatal_message("Couldn't create detect variants dispatcher from params:\n".Data::Dumper::Dumper \%params);
    }
    my $rv = $command->execute;
    unless ($rv){
        $self->fatal_message("Failed to execute detect variants dispatcher with params:\n".Data::Dumper::Dumper \%params);
    }

    $self->debug_message("detect variants command completed successfully");

    return 1;
}

sub _should_skip_run {
    my $self = shift;
    my $pp = $self->build->processing_profile;

    if(defined $pp->snv_detection_strategy || defined $pp->indel_detection_strategy ||
            defined $pp->sv_detection_strategy || defined $pp->cnv_detection_strategy) {
        return;
    }
    else {
        $self->debug_message('No strategies defined--decided to skip detect variants step.');
        return 1;
    }
}

1;
