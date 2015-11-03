package Genome::Model::SingleSampleGenotype::Command::PrepareReferenceBuckets;

use strict;
use warnings;

use Genome;

class Genome::Model::SingleSampleGenotype::Command::PrepareReferenceBuckets {
    is => 'Genome::Model::SingleSampleGenotype::Command::Base',
    doc => 'Prepare buckets for running the HaplotypeCaller in parallel',
};

sub _result_accessor {
    return 'output_result';
}

sub _command_class {
    return 'Genome::Model::ReferenceSequence::Command::CreateBuckets';
}

sub _params_for_command {
    my $self = shift;
    my $build = $self->build;

    my $result_users = Genome::SoftwareResult::User->user_hash_for_build($build);

    my @params = (
        reference_sequence_build => $build->reference_sequence_build,
        result_users => $result_users,
    );

    return @params;
}

sub sub_command_category {
    return 'analyst tools';
}

1;
