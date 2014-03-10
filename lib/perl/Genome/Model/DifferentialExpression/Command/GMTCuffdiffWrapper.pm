# FIXME This is only here until unittests around Genome::Command::WithSoftwareResult
# are in place.
package Genome::Model::DifferentialExpression::Command::GMTCuffdiffWrapper;

use strict;
use warnings;

use Genome;
use Cwd 'abs_path';


class Genome::Model::DifferentialExpression::Command::GMTCuffdiffWrapper {
    is => ['Genome::Command::WithSoftwareResult'],
    has => [
        _software_result_version => { # required by WithSoftwareResult
            is => 'Integer',
            is_param => 1,
            valid_values => [1],
            default_value => '1',
            doc => 'the version of results, which may iterate as execute_logic iterates',
        },
    ],
    has_input => [
        transcript_gtf_file => {
            doc => 'A GTF format file of transcript annotation to perform differential expression tests.',
        },
        condition_model_ids_string => {
            is => 'String',
            doc => 'A list of (model_id lists [comma separated]) that are space separated',
        },
    ],
    has_transient => [
        output_directory => {
            doc => 'The directory to write output files.',
        },
    ],
    has_param => [
        use_version => {
            is => 'Version',
            is_optional => 1,
        },
        cuffdiff_params => {
            is_optional => 1,
            doc => 'A text string of optional parameters to pass to cuffdiff',
        },
    ],
};

sub software_result_type {
    return 'staged';
}

sub _execute {
    my ($self, $output_dir) = @_;

    $self->debug_message("Executing cuffdiff with output_directory set to: $output_dir\n");
    if ($self->_execute_gmt_cuffdiff($output_dir)) {
        return 1;
    } else {
        $self->error_message("Failed to run Genome::Model::Tools::Cufflinks::Cuffdiff->execute");
        return 0;
    }
}

# Called after result _promote_data()
sub _finalize {
    my ($self, $result) = @_;

    my $output_directory = $result->output_dir;
    $self->debug_message(
        sprintf("Symlinking results from %s to %s",
            $output_directory, $self->output_directory)
    );
    Genome::Sys->symlink_directory($output_directory, $self->output_directory);
}

sub _execute_gmt_cuffdiff {
    my $self = shift;
    my $output_directory = shift;

    my $bam_file_paths = _resolve_bam_file_paths($self->condition_model_ids_string);

    my $cmd = Genome::Model::Tools::Cufflinks::Cuffdiff->create(
        params => $self->cuffdiff_params,
        output_directory => $output_directory,
        bam_file_paths => $bam_file_paths,
        transcript_gtf_file => $self->transcript_gtf_file,
        use_version => $self->use_version,
    );
    return $cmd->execute();
}

sub _resolve_bam_file_paths {
    my ($condition_model_ids_string) = @_;

    my @condition_model_ids = split(/ /, $condition_model_ids_string);
    my $bam_file_paths = '';
    for my $condition_model_ids (@condition_model_ids) {
        my @model_ids = split(/,/,$condition_model_ids);
        my @condition_bam_file_paths;
        for my $model_id (@model_ids) {
            my $bam_file = _get_bam_file($model_id);
            push @condition_bam_file_paths, $bam_file;
        }
        my $condition_bam_files_string = join(',',@condition_bam_file_paths);
        $bam_file_paths .= $condition_bam_files_string .' ';
    }
    return $bam_file_paths;
}

sub _get_bam_file {
    my $model_id = shift;

    my $rna_seq_model = Genome::Model->get($model_id);
    # Or should we just fine the latest AlignmentResult....
    my $last_succeeded_build = $rna_seq_model->last_succeeded_build;
    my $bam_file = abs_path($last_succeeded_build->alignment_result->bam_file);

    return $bam_file;
}

1;
