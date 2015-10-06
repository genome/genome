package Genome::Model::ReferenceVariation::Command::HaplotypeCaller;

use strict;
use warnings;

use Genome;

class Genome::Model::ReferenceVariation::Command::HaplotypeCaller {
    is => 'Genome::Model::ReferenceVariation::Command::Base',
    doc => 'Runs the Haplotype Caller on the alignments',
    has_optional_input => [
        intervals => {
            is => 'Text',
            doc => 'intervals to restrict this run of the haplotype caller',
        },
    ],
};

sub _result_accessor {
    return 'output_result';
}

sub _command_class {
    return 'Genome::Model::ReferenceVariation::Result::HaplotypeCallerWrapper';
}

sub _params_for_command {
    my $self = shift;
    my $build = $self->build;

    my $result_users = Genome::SoftwareResult::User->user_hash_for_build($build);
    $result_users->{haplotype_caller_result} = $build;

    my @params = (
        alignment_result => $build->merged_alignment_result,
        result_users => $result_users,
        haplotype_caller_version => $build->model->haplotype_caller_version,
        emit_reference_confidence => 'GVCF',
        gvcf_gq_bands => [5,20,60],
        intervals => $self->intervals,
    );

    return @params;
}

1;
