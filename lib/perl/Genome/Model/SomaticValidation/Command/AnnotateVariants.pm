package Genome::Model::SomaticValidation::Command::AnnotateVariants;

use strict;
use warnings;

use Genome;

use Data::Dumper qw();

use File::Copy qw();
use File::Spec qw();
use File::Basename qw();

class Genome::Model::SomaticValidation::Command::AnnotateVariants {
    is => 'Genome::Command::Base',
    has => [
        annotator_version => {
            doc => 'Version of the "annotate transcript-variants" tool to run   during the annotation step',
            default_value => Genome::Model::Tools::Annotate::TranscriptVariants->default_annotator_version,
            is_optional => 1,
            is_input => 1,
        },
        build => {
            is => 'Genome::Model::Build::SomaticValidation',
            id_by => 'build_id',
        },
        build_id => {
            is => 'Integer',
            is_input => 1,
            is_output => 1,
            doc => 'build id of SomaticValidation model',
        },
    ],
    has_param => [
        lsf_queue => {
            default => 'apipe',
        },
        snv_tiers_to_annotate => {
            is => 'Array',
            is_input => 1,
            is_optional => 1,
            default_value => [1],
        },
    ],
};

sub execute {
    my $self = shift;
    my $build = $self->build;

    $self->status_message("Executing Annotate step");


    my $annotation_directory = $build->data_directory . '/annotations';
    $self->status_message(sprintf("Creating annotation directory at %s",
            $annotation_directory));
    Genome::Sys->create_directory($annotation_directory);

    $self->status_message("Annotating SNVs");
    $self->annotate_snvs($annotation_directory);

    $self->status_message("Annotation complete.");

    return 1;
}

sub annotate_snvs {
    my ($self, $annotation_directory) = @_;

    for my $tier (@{$self->snv_tiers_to_annotate}) {
        TYPE:
        for my $type (qw(hq lq)) {
            $self->status_message(sprintf("Annotating tier %s %s SNVs", $tier, $type));

            my $temp_output_path = Genome::Sys->create_temp_file_path();

            my $tiering_result = $self->build->result_user(
                label => sprintf('snv_tiered_%s', $type))->software_result;
            my $input_path = $tiering_result->path(sprintf('*tier%s*.bed', $tier));
            $self->status_message(sprintf('Preparing to annotate file: %s',
                    $input_path));

            unless (defined $input_path) {
                die $self->error_message(sprintf(
                        'Could not resolve path for %s input bed file for tier %s.',
                        $type, $tier));
            }
            unless (-e $input_path) {
                $self->error_message(sprintf(
                        'Skipping annotation for non-existant input file %s',
                        $input_path));
                next TYPE;
            }

            my $output_path = $self->_resolve_output_path($annotation_directory,
                                                          $input_path);

            unless (-s $input_path) {
                $self->status_message(sprintf('Skipping annotation for empty file %s',
                        $input_path));
                Genome::Sys->write_file($output_path, '');
                next TYPE;
            }

            my %annotation_params = (
                no_headers => 1,
                use_version => $self->annotator_version,
                build_id => $self->build->annotation_build->id,
                variant_bed_file => $input_path,
                annotation_filter => 'top',
                output_file => $temp_output_path,
            );

            $self->status_message(sprintf(
                    'Executing GMT:AnnotateTranscriptVariants with parameters: %s',
                    join('', Data::Dumper::Dumper(%annotation_params))));

            my $annotation_rv = Genome::Model::Tools::Annotate::TranscriptVariants->execute(
                %annotation_params);
            my $annotation_error = $@;

            unless ($annotation_rv) {
                die $self->error_message(sprintf(
                        'Failed to execute GMT:Annotate::TranscriptVariants: %s', $@));
            }

            # verify output file exists?
            unless (-s $temp_output_path) {
                die $self->error_message(sprintf(
                        'No output file found from GMT:Annotate::TranscriptVariants for tier %s SNVs.',
                        $tier));
            }

            File::Copy::copy($temp_output_path, $output_path);
        }
    }
}

sub _resolve_output_path {
    my ($self, $annotation_directory, $input_path) = @_;

    my $filename = File::Basename::basename($input_path);

    return File::Spec->join($annotation_directory, $filename);
}

1;
