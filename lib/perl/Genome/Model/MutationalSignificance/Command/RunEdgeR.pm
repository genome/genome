package Genome::Model::MutationalSignificance::Command::RunEdgeR;

use strict;
use warnings;
use Genome;

class Genome::Model::MutationalSignificance::Command::RunEdgeR {
    is => 'Command::V2',
    has_input => [
        clinseq_models => {
            is => 'Genome::Model::ClinSeq',
            is_many => 1,
            is_input => 1,
        },
        counts_per_million => {
            is => 'Integer',
        },
        num_samples => {
            is => 'Integer',
        },
        p_value => {
            is => "Number",
            doc => "The p-value used for significance testing (0 < p < 1)",
            default_value => 0.05,
        },
    ],
    has_output => [
        output_file => {
            is => "Text",
            doc => "Output file path",
        },
    ],
};

sub execute {
    my $self = shift;

    my $counts_file = Genome::Sys->create_temp_file_path();

    my $make_gene_counts_file = Genome::Model::MutationalSignificance::Command::MakeGeneCountsFile->execute(
        clinseq_models  => [$self->clinseq_models],
        counts_file     => $counts_file,
    );

    my $filtered_counts_file = Genome::Sys->create_temp_file_path();

    my $filter_counts_file = Genome::Model::Tools::EdgeR::FilterCountsFile->execute(
        counts_file             => $counts_file,
        counts_per_million      => $self->counts_per_million,
        num_samples             => $self->num_samples,
        filtered_counts_file    => $filtered_counts_file,
    );

    my $edge_r = Genome::Model::Tools::EdgeR::Glm->execute(
        counts_file => $filter_counts_file->filtered_counts_file,
        groups      => $make_gene_counts_file->groups,
        subjects    => $make_gene_counts_file->subjects,
        output_file => $self->output_file,
        p_value     => $self->p_value,
    );
}

1;
