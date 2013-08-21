package Genome::Model::SomaticValidation::Command::ProcessValidation;

use strict;
use warnings;

use Genome;
use Cwd qw(abs_path);

class Genome::Model::SomaticValidation::Command::ProcessValidation {
    is => 'Command',
    has_input => [
        filtered_validation_file    => { is => 'Text', doc => "bed file of variants passing filter", is_optional => 0 },
        min_coverage                => { is => 'Text', doc => "Minimum coverage to call a site", is_optional => 1 },
        output_file                 => { is => 'Text', doc => "Output file for validation results", is_output => 1 },
        output_plot                 => { is => 'Boolean', doc => "Optional plot of variant allele frequencies", is_optional => 1, },
        build_id => {
            is => 'Integer',
            is_output => 1,
            doc => 'build id of SomaticValidation model',
        },
    ],
    has => [
        build => {
            is => 'Genome::Model::Build::SomaticValidation',
            id_by => 'build_id',
        },
    ],
    has_param => [
        lsf_resource => {
            default_value => "-M 30000000 -R 'rusage[tmp=2000 && mem=30000] select[mem>30000 && tmp>2000]'",
        },
    ],
    doc => 'final processing of HQ SNV detection results',
};

sub sub_command_category { 'pipeline steps' }

sub shortcut {
    my $self = shift;
    my $build = $self->build;

    my $variant_list = $build->snv_variant_list;
    unless($variant_list) {
        $self->status_message('No SNVs list provided. Skipping run.');
        return 1;
    }

    return; #have to do the work otherwise
}

sub execute {
    my $self = shift;
    my $build = $self->build;

    my $variant_list = $build->snv_variant_list;
    unless($variant_list) {
        $self->status_message('No SNVs list provided. Skipping run.');
        return 1;
    }

    my ($snv_variant_file) = glob($variant_list->output_dir . '/snvs.hq.bed');
    unless($snv_variant_file) {
        $self->error_message('Failed to get a snv variant file for this build!');
        die $self->error_message;
    }

    my $anno_file = $build->data_directory . '/variant_list.snvs';
    my $bed_to_anno = Genome::Model::Tools::Bed::Convert::BedToAnnotation->create(
        snv_file => $snv_variant_file,
        output => $anno_file,
    );
    unless($bed_to_anno->execute) {
        die $self->error_message('Failed to convert BED file to annotation format');
    }

    my @validation_original_file = glob($build->data_directory . '/variants/snv/varscan-somatic-validation*/snvs.hq.validation'); 
    unless(scalar @validation_original_file == 1) {
        die $self->error_message('Unable to determine the original varscan file to use for ProcessValidation run');
    }

    my $filtered_bed = $self->filtered_validation_file;
    my $filtered_original_file = Cwd::abs_path($filtered_bed);
    $filtered_original_file =~ s/(?:\.v\d+)?\.bed$//;
    unless(Genome::Sys->check_for_path_existence($filtered_original_file)) {
        $self->error_message('Failed to find original filtered file to use for ProcessValidation run. This step requires that the final SNV result contain a Varscan output.');
        return 1;
    }

    my $process_validation = Genome::Model::Tools::Varscan::ProcessValidation->create(
        validation_file => $validation_original_file[0],
        filtered_validation_file => $filtered_original_file,
        variants_file => $anno_file,
        output_file => $self->output_file,
        output_plot => $self->output_plot,
    );

    unless($process_validation->execute) {
        die $self->error_message('Execution of ProcessValidation failed');
    }

    return 1;
}

1;
