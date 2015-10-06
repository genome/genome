package Genome::Model::ReferenceVariation::Command::AlignReads;

use strict;
use warnings;

use Genome;

class Genome::Model::ReferenceVariation::Command::AlignReads {
    is => 'Command::V2',
    has_input => [
        build => {
            is => 'Genome::Model::Build::ReferenceVariation',
            id_by => 'build_id',
            doc => 'build for which to run alignments',
            is_output => 1,
        },
    ],
    doc => 'Align the instrument data to the reference for the build.',
};

sub shortcut {
    my $self = shift;

    my $cmd = $self->_alignment_command;
    my $retval = $cmd->shortcut;

    if ($retval and $cmd->result_id) {
        $self->status_message('Found existing result: %s', $cmd->result_id);
    }

    return $retval;
}

sub execute {
    my $self = shift;

    my $cmd = $self->_alignment_command;
    my $retval = $cmd->execute;

    if ($retval and $cmd->result_id) {
        $self->status_message('Generated result: %s', $cmd->result_id);
    } else {
        die $self->error_message('Failed to produce aligned BAM.');
    }

    return $retval;
}

sub _alignment_command {
    my $self = shift;

    my @params = $self->_params_for_alignment;

    my $cmd = Genome::InstrumentData::Command::AlignAndMerge->create(
        @params,
    );

    die $self->error_messaage('Failed to create alignment command!') unless $cmd;

    return $cmd;
}

sub _params_for_alignment {
    my $self = shift;
    my $build = $self->build;

    my $result_users = Genome::SoftwareResult::User->user_hash_for_build($build);
    $result_users->{merged_alignment_result} = $build;

    my $api_inputs = Genome::InstrumentData::Composite::Workflow::Generator->inputs_for_api_version(
        $build->aligner_api_version
    );

    my @params = (
        instrument_data => [$build->instrument_data],
        reference_sequence_build => $build->reference_sequence_build,
        result_users => $result_users,
        name => 'speedseq',
        version => $build->aligner_version,
        params => $build->aligner_params,
        %$api_inputs,
    );

    return @params;
}

1;
