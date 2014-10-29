package Genome::Model::DifferentialExpression::Command::Cummerbund;

use strict;
use warnings;

use Genome;

my $DEFAULT_LSF_RESOURCE = "-R 'select[mem>=16000] rusage[mem=16000] span[hosts=1]' -M 16000000 ";

class Genome::Model::DifferentialExpression::Command::Cummerbund {
    is => ['Command::V2'],
    has => [
        build => { is => 'Genome::Model::Build', id_by => 'build_id', },
    ],
    has_input_output => {
        build_id => {},
    },
    has_param => [
        lsf_resource => { default_value => $DEFAULT_LSF_RESOURCE },
    ],
};

sub execute {
    my $self = shift;

    my $output_dir = $self->build->summarize_differential_expression_directory;
    unless (-d $output_dir) {
        Genome::Sys->create_directory($output_dir);
    }
    my $cmd = Genome::Model::Tools::Cufflinks::Cummerbund->execute(
        cuffdiff_directory => $self->build->differential_expression_directory,
        pdf_report_file => $self->build->summary_report_pdf_file,
    );
    unless ($cmd and $cmd->result) {
        $self->error_message('Failed to run cummerbund!');
        return;
    }
    return 1;
}

1;

