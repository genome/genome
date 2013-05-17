package Genome::Model::DifferentialExpression::Command::GMTCuffdiffWrapper;

use strict;
use warnings;

use Genome;


class Genome::Model::DifferentialExpression::Command::GMTCuffdiffWrapper {
    is => ['Genome::Command::WithSavedResults'],
    has_input => [
        transcript_gtf_file => {
            doc => 'A GTF format file of transcript annotation to perform differential expression tests.',
        },
        bam_file_paths => {
            doc => 'A string of bam file paths separated by white space. Supply replicate SAMs as comma separated lists for each condition: sample1_rep1.sam,sample1_rep2.sam,...sample1_repM.sam',
        },
    ],
    has => [
        output_directory => {
            doc => 'The directory to write output files.',
        },
        stage_output => {
            default_value => 1,
            doc => 'Should this command stage its result in temp until it has completed?',
        },
    ],
    has_param => [
        use_version => {
            is => 'Version',
            is_optional => 1,
        },
        result_version => {
            is => 'Integer',
            valid_values => [1],
            default_value => '1',
            doc => 'the version of results, which may iterate as execute_logic iterates',
        },
        cuffdiff_params => {
            is_optional => 1,
            doc => 'A text string of optional parameters to pass to cuffdiff',
        },
    ],
};

sub _execute_v1 {
    my $self = shift;

    # $self->output_dir is set up by WithSavedResults and value set to ->result->temp_staging_directory
    # because ->stage_output = 1 and ->output_dir = undef before _execute_v1 is called.
    # I don't like that this comment is needed (means things are not obvious).
    my $output_directory = $self->output_dir;
    $self->status_message("Setting cuffdiff's output_directory to: $output_directory\n");

    if ($self->_execute_gmt_cuffdiff($output_directory)) {
        $self->status_message(
            sprintf("Symlinking results from %s to %s",
                $output_directory, $self->output_directory)
        );
        Genome::Sys->symlink_directory($output_directory, $self->output_directory);
        return 1;
    } else {
        $self->error_message("Failed to run Genome::Model::Tools::Cufflinks::Cuffdiff->execute");
        return 0;
    }
}

sub _execute_gmt_cuffdiff {
    my $self = shift;
    my $output_directory = shift;

    my $cmd = Genome::Model::Tools::Cufflinks::Cuffdiff->get(
        params => $self->params,
        output_directory => $output_directory,
        bam_file_paths => $self->bam_file_paths,
        transcript_gtf_file => $self->transcript_gtf_file,
        use_version => $self->use_version,
    );
    return $cmd->execute();
}
