package Genome::Model::ReferenceVariation::Command::AlignReads;

use strict;
use warnings;

use Genome;

class Genome::Model::ReferenceVariation::Command::AlignReads {
    is => 'Genome::Model::ReferenceVariation::Command::Base',
    doc => 'Align the instrument data to the reference for the build.',
};

sub _result_accessor {
    return 'alignment_result';
}

sub _command_class {
    return 'Genome::InstrumentData::Command::AlignAndMerge';
}

sub _params_for_command {
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
        name => $build->aligner_name,
        version => $build->aligner_version,
        params => $build->aligner_params,
        %$api_inputs,
    );

    return @params;
}

1;
