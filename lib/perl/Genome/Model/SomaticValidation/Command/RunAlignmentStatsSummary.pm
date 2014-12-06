package Genome::Model::SomaticValidation::Command::RunAlignmentStatsSummary;

use strict;
use warnings;

use Genome;

class Genome::Model::SomaticValidation::Command::RunAlignmentStatsSummary {
    is => 'Genome::Command::DelegatesToResult',
    doc => 'Generate a spreadsheet, tsv file, of alignment metrics for Somatic Validation models.  Duplicate samples using the same processing profile will only be reported once.',
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
    has_optional_input => [
        haploid_coverage => {
            is => 'Boolean',
            doc => 'Execute bam-check to calculate the haploid coverage. This can take a LONG time!',
            default_value => 0,
        },
        targeted_insert_length => {
            is => 'Boolean',
            doc => 'Add the targeted insert length for each sequence item.',
            default_value => 0,
        },
    ],
};

sub help_detail {
    return "Summarize the alignment stats for all listed somatic validation models.  One line will be output in the tsv file for each sample."
}

sub result_class {
    return "Genome::Model::SomaticValidation::Command::AlignmentStatsSummary";
}

sub input_hash {
    my $self = shift;

    my %hash = (
        builds => [$self->_builds],
        haploid_coverage => $self->haploid_coverage,
        targeted_insert_length => $self->targeted_insert_length,
    );

    return %hash;
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
