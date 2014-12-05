package Genome::Model::SomaticValidation::Command::RunCoverageStatsSummary;

use strict;
use warnings;

use Genome;

class Genome::Model::SomaticValidation::Command::RunCoverageStatsSummary {
    is => 'Genome::Command::DelegatesToResult',
    doc => 'Generate a spreadsheet, tsv file, of coverage metrics for Somatic Validation models.  Duplicate normal samples will only be reported once.',
    has_input => [
        models => {
            is => 'Genome::Model::SomaticValidation',
            is_many => 1,
            shell_args_position => 1,
            doc => 'The Somatic Validation models or an expression to resolve the Somatic Validation models.',
        },
        user => {
            is => 'Genome::Process',
            is_optional => 1,
            id_by => 'process_id',
        },
        process_id => {
            is => 'Text',
        },
    ],
};

sub help_detail {
    return "Summarize the coverage stats for all listed somatic validation builds.  One line will be output in the tsv file for each sample, region of interest, wingspan and minimum depth filter."
}

sub result_class {
    return "Genome::Model::SomaticValidation::Command::CoverageStatsSummary";
}

sub input_hash {
    my $self = shift;

    return (
        builds => [$self->_builds],
    );
}

sub _builds {
    my $self = shift;
    my @builds;
    for my $model ($self->models) {
        unless (defined $model->last_succeeded_build) {
            die sprintf("No succeeded build found for model %s", $model->id);
        }
        push @builds, $model->last_succeeded_build;
    }
    return @builds;
}
1;
