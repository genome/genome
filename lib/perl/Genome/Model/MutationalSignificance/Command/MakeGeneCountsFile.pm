package Genome::Model::MutationalSignificance::Command::MakeGeneCountsFile;

use strict;
use warnings;
use Genome;

class Genome::Model::MutationalSignificance::Command::MakeGeneCountsFile {
    is => 'Command::V2',
    has_input => [
        clinseq_models => {
            is => 'Genome::Model::ClinSeq',
            is_many => 1,
            is_input => 1,
        },
    ],
    has_output => [
        counts_file => {
            is => 'File',
        },
        subjects => {
            is => 'Text',
            is_optional => 1,
        },
        groups => {
            is => 'Text',
            is_optional => 1,
        }
    ],
};

sub execute {
    my $self = shift;

    my @clinseq_models = $self->clinseq_models;

    my (@tumor_models, @normal_models);

    for my $model (@clinseq_models) {
        if ( defined($model->tumor_rnaseq_model)) {
                push @tumor_models, $model->tumor_rnaseq_model;
        }
        if ( defined($model->normal_rnaseq_model) ) {
                push @normal_models, $model->normal_rnaseq_model;
        }
    }

    my $cmd = Genome::Model::RnaSeq::Command::MakeGeneCountsFile->create(
        normal_models => \@normal_models,
        tumor_models => \@tumor_models,
        counts_file => $self->counts_file,
    );

    $cmd->execute;

    $self->counts_file($cmd->counts_file);
    $self->subjects($cmd->subjects);
    $self->groups($cmd->groups);

    return 1;
}

1;
