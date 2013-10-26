package Genome::Model::SomaticValidation::Command::ExtractVariantsToValidate;

use strict;
use warnings;

use Genome;
use Genome::Model::Tools::DetectVariants2::Utilities qw(
    final_result_for_variant_type
);

class Genome::Model::SomaticValidation::Command::ExtractVariantsToValidate {
    is => 'Command::V2',
    has_input => [
        model => {
            is => 'Genome::Model',
            doc => 'the somatic model from discovery',
        },
        target => {
            is => 'Genome::FeatureList',
            doc => 'The target region set being used for validation',
        },
        report_only => {
            is => 'Boolean',
            is_optional => 1,
            default_value => 0,
            doc => 'if specified will not create the manual results but will only report statistics about the coverage of the previously detected variants',
        },
    ],
    has_optional_output => {
        snv_variant_list => {
            is => 'Genome::Model::Tools::DetectVariants2::Result::Manual',
            doc => 'the list of SNVs from discovery intersected with the target set',
        },
        indel_variant_list => {
            is => 'Genome::Model::Tools::DetectVariants2::Result::Manual',
            doc => 'the list of indels from discovery intersected with the target set',
        },
        sv_variant_list => {
            is => 'Genome::Model::Tools::DetectVariants2::Result::Manual',
            doc => 'the list of SVs from discovery intersected with the target set',
        },
        cnv_variant_list => {
            is => 'Genome::Model::Tools::DetectVariants2::Result::Manual',
            doc => 'the list of CNVs from discovery intersected with the target set',
        },
    },
};

sub sub_command_category { 'analyst tools' }

sub execute {
    my $self = shift;

    my $build = $self->model->last_complete_build;
    unless($build) {
        die $self->error_message('No complete build found on the model.');
    }

    my $target = $self->target;

    VARIANT_TYPE: for my $variant_type ('snv', 'indel') { #('snv', 'indel', 'sv', 'cnv') { #TODO support other variant types?

        my @build_results = $build->results;
        my @results = grep($_->isa('Genome::Model::Tools::DetectVariants2::Classify::Tier') && $_->variant_type eq $variant_type, @build_results);

        if(!@results) {
            @results = final_result_for_variant_type(\@build_results, $variant_type);
        }

        unless(@results) {
            $self->warning_message('No result found for ' . $variant_type . ' to extract validation variants');
            next VARIANT_TYPE;
        }

        my @output_files;
        for my $result (@results) {
            my $output_dir = $result->output_dir;
            my @f = glob($output_dir . '/*.bed');

            for my $f (@f) {
                my ($hit, $miss) = $self->_intersect_file_with_targets($f);
                push @output_files, $hit;
                my $hit_count = `wc -l $hit`;
                ($hit_count) = $hit_count =~ m/^(\d+)/;
                my $miss_count = `wc -l $miss`;
                ($miss_count) = $miss_count =~ m/^(\d+)/;

                $self->status_message('Found ' . $hit_count . ' variants in the target set for ' . $f . '. (' . sprintf("%.3f", ($hit_count / ($hit_count+$miss_count))*100) . '%)');
            }
        }

        unless($self->report_only) {
            my $final_file = Genome::Sys->create_temp_file_path;
            my $merge_cmd = Genome::Model::Tools::Joinx::Sort->create(
                input_files => \@output_files,
                merge_only => 1,
                output_file => $final_file,
            );

            unless($merge_cmd->execute()) {
                die $self->error_message('Failed to generate merged file.');
            }

            my $manual_result_cmd = Genome::Model::SomaticValidation::Command::ManualResult->create(
                source_build => $build,
                variant_file => $final_file,
                variant_type => $variant_type,
                description => 'subset of variants from build ' . $build->id . ' contained within the target set ' . $target->name,
                format => 'bed',
            );

            unless($manual_result_cmd->execute()) {
                die $self->error_message('Failed to generate manual result');
            }

            $self->status_message('Generated manual result for ' . $variant_type . ': ' . $manual_result_cmd->manual_result->id);
            my $accessor = $variant_type . '_variant_list';
            $self->$accessor($manual_result_cmd->manual_result);
        }
    }

    return 1;
}

sub _intersect_file_with_targets {
    my $self = shift;
    my $f = shift;
    my $target = $self->target;
    my $target_file = $target->file_path;

    my $miss_a = Genome::Sys->create_temp_file_path;
    my $output = Genome::Sys->create_temp_file_path;

    my $intersect_cmd = Genome::Model::Tools::Joinx::Intersect->create(
        input_file_a => $f,
        input_file_b => $target_file,
        miss_a_file => $miss_a,
        output_file => $output,
    );

    unless($intersect_cmd->execute()) {
        die $self->error_message('Failed to execute joinx intersection using ' . $f . ' and ' . $target_file);
    }

    return ($output, $miss_a);
}

1;

