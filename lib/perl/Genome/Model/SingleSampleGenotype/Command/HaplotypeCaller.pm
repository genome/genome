package Genome::Model::SingleSampleGenotype::Command::HaplotypeCaller;

use strict;
use warnings;

use Genome;

class Genome::Model::SingleSampleGenotype::Command::HaplotypeCaller {
    is => 'Genome::Model::SingleSampleGenotype::Command::Base',
    doc => 'Runs the Haplotype Caller on the alignments',
    has_optional_input => [
        intervals => {
            is => 'Text',
            doc => 'intervals to restrict this run of the haplotype caller',
            is_many => 1,
        },
    ],
};

sub _result_accessor {
    return 'output_result';
}

sub _command_class {
    return 'Genome::Model::SingleSampleGenotype::Result::HaplotypeCallerWrapper';
}

sub _params_for_command {
    my $self = shift;
    my $build = $self->build;

    my $result_users = Genome::SoftwareResult::User->user_hash_for_build($build);

    my @params = (
        alignment_result => $build->merged_alignment_result,
        result_users => $result_users,
        haplotype_caller_version => $build->model->haplotype_caller_version,
        emit_reference_confidence => 'GVCF',
        gvcf_gq_bands => [5,20,60],
    );

    my $intervals = $self->intervals;
    push @params, intervals => $intervals if $intervals;

    return @params;
}

sub _label_for_result {
    return 'haplotype_caller_result';
}

1;
