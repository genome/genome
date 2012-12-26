package Genome::Model::SomaticValidation::Command::ValidateLargeIndels;

use strict;
use warnings;

use Genome;

use Data::UUID qw();
use File::Spec qw();

class Genome::Model::SomaticValidation::Command::ValidateLargeIndels {
    is => 'Genome::Command::Base',
    has => [
        build => {
            is => 'Genome::Model::Build::SomaticValidation',
            id_by => 'build_id',
        },
        build_id => {
            is => 'Integer',
            is_input => 1,
            doc => 'build id of SomaticValidation model',
        },
        _long_indel_bed_file => {
            is => 'Text',
            is_input => 1,
            is_optional => 1,
            doc => 'Path to bed file that contains the large indels. Resolved via the build.',
        },
        reference_transcripts => {
            is => 'Text',
            default => "NCBI-human.ensembl/67_37l_v2", #TODO this should be a param from the somatic validation processing profile
            doc => 'The set of reference_transcripts to use, which get passed to the annotator',
        },
    ],
    has_param => [
        lsf_queue => {
            default => 'apipe',
        },
    ],
};

sub sub_command_category { 'pipeline steps' }

sub execute {
    my $self = shift;

    return 1 unless $self->build->normal_sample;

    my $long_indel_bed_file = $self->_resolve_long_indel_bed_file;
    unless ($long_indel_bed_file) {
        $self->warning_message("No long indel bed file exists with size. Skipping validation.");
        return 1;
    }

    $self->_long_indel_bed_file($long_indel_bed_file);
    my $output_directory = $self->_create_output_directory();

    $self->status_message("Running step 1");
    my ($tumor_model, $normal_model) = $self->_run_indel_step_1(
            $output_directory);

    $self->status_message("Starting and waiting for builds");
    $self->_start_and_wait_for_builds($tumor_model, $normal_model);

    $self->status_message("Running step 2");
    $self->_run_indel_step_2($output_directory, $tumor_model, $normal_model);

    return 1;
}

sub _resolve_output_directory {
    my $self = shift;
    return File::Spec->join($self->build->data_directory, 'large_indel_validation');
}

sub _resolve_long_indel_bed_file {
    my $self = shift;
    my $long_indel_bed_file = $self->build->data_directory . "/indel_validation/large_indels.bed";
    unless (-s $long_indel_bed_file) {
        $self->warning_message("Long indel bed file $long_indel_bed_file does not exist or has no size");
        return;
    }
    return $long_indel_bed_file;
}

sub _create_output_directory {
    my $self = shift;

    my $output_directory = $self->_resolve_output_directory();
    Genome::Sys->create_directory($output_directory);

    return $output_directory;
}

sub _run_indel_step_1 {
    my ($self, $output_directory) = @_;

    my $sample_identifier = Data::UUID->new->create_str();
    #TODO the instructions say to "be sure to save the STDOUT from this tool

    my $command = Genome::Model::Tools::Validation::LongIndelsPartOne->create(
        somatic_validation_model_id => $self->build->model->id,
        long_indel_bed_file => $self->_long_indel_bed_file,
        sample_identifier => $sample_identifier,
        output_dir => $output_directory,
        reference_transcripts => $self->reference_transcripts,
    );

    unless ($command->execute) {
        die $self->error_message("Failed to execute LongIndelsPartOne");
    }

    return ($command->_new_tumor_model, $command->_new_normal_model);
}

sub _start_and_wait_for_builds {
    my ($self, @models) = @_;
    local $SIG{CHLD} = 'DEFAULT'; # FIXME this should not have to be here, really.
    my $command = Genome::Model::Build::Command::Start->create(models => \@models);
    unless ($command->execute) {
        die $self->error_message("Failed to start builds");
    }
    UR::Context->commit;

    my @alignment_builds = $command->builds;
    my @alignment_build_ids = map{$_->id} @alignment_builds;
    $self->_wait_for_builds(@alignment_build_ids);

    return @alignment_build_ids;
}

sub _wait_for_builds {
    my ($self, @alignment_build_ids) = @_;

    my $max_sleeps = 288; # Wait for 48 hours max;
    my $current_sleeps = 0;
    for my $build_id (@alignment_build_ids) {
        my $build = UR::Context->current->reload("Genome::Model::Build", id => $build_id);
        UR::Context->current->reload($build->the_master_event);

        if ($build->status eq "Running" || $build->status eq 'Scheduled') {
            $self->status_message("Build " . $build->id . " has a status of " . $build->status . " ... waiting 10 minutes for it to finish.");
            $current_sleeps++;
            if ($current_sleeps > $max_sleeps) {
                die $self->error_messsage("Slept for 48 hours. THAT FELT GREAT. But we've probably been waiting for these builds for too long so I'm dying");
            }
            sleep 600;
            redo;
        }
        unless ($build->status eq "Succeeded") {
            die $self->error_message("Build " . $build->id . " has a status of " . $build->status . " ... cannot continue.");
        }
    }

    return 1;
}

sub _run_indel_step_2 {
    my ($self, $output_directory, $tumor_model, $normal_model) = @_;

    my $annotation_build = Genome::Model::Build::ImportedAnnotation->get(name => $self->reference_transcripts);
    
    return Genome::Model::Tools::Validation::LongIndelsPartTwo->execute(
        output_dir => $output_directory,
        tumor_val_model_copy_id => $tumor_model->id,
        normal_val_model_copy_id => $normal_model->id,
        contigs_file => "$output_directory/contigs.fa",
        tier_file_location => $annotation_build->tiering_bed_files_by_version(3), 
    );
}

1;
