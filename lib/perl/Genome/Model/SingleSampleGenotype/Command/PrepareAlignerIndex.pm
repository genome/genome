package Genome::Model::SingleSampleGenotype::Command::PrepareAlignerIndex;

use strict;
use warnings;

use Genome;

class Genome::Model::SingleSampleGenotype::Command::PrepareAlignerIndex {
    is => 'Genome::Model::SingleSampleGenotype::Command::Base',
    doc => 'Prepare an aligner index for use in the single-sample genotype pipeline',
    has_param => [
        lsf_resource => {
            is => 'Text',
            value => __PACKAGE__->_command_class->__meta__->property(property_name => 'lsf_resource')->default_value,
        },
    ],
};

sub _result_accessor {
    return 'aligner_index';
}

sub _command_class {
    return 'Genome::Model::ReferenceSequence::Command::CreateAlignerIndex';
}

sub _params_for_command {
    my $self = shift;
    my $build = $self->build;
    my $model = $build->model;

    my $result_users = Genome::SoftwareResult::User->user_hash_for_build($build);

    my @params = (
        aligner_name => $model->aligner_name,
        aligner_params => $model->aligner_params || '',
        aligner_version => $model->aligner_version,
        reference_sequence_build => $build->reference_sequence_build,
        result_users => $result_users,
    );

    return @params;
}

sub sub_command_category {
    return 'analyst tools';
}

1;
