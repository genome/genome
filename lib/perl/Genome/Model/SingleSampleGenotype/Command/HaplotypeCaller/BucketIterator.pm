package Genome::Model::SingleSampleGenotype::Command::HaplotypeCaller::BucketIterator;

use strict;
use warnings;

use Genome;

class Genome::Model::SingleSampleGenotype::Command::HaplotypeCaller::BucketIterator {
    is => 'Command::V2',
    doc => 'Runs the Haplotype Caller on the alignments',
    has_input => [
        bucket => {
            is => 'Text',
            doc => 'bucket to restrict this run of the haplotype caller',
        },
        build => {
            is => 'Genome::Model::Build::SingleSampleGenotype',
            doc => 'build for which to run the command',
            is_output => 1,
            shell_args_position => 1,
        },
    ],
    has_param => [
        lsf_queue => {
            is => 'Text',
            default => Genome::Config::get('lsf_queue_build_worker_alt'),
        },
    ],
    doc => 'Runs the haplotype-caller command for each interval in the bucket',
};

sub shortcut {
    return $_[0]->_run('shortcut');
}

sub execute {
    return $_[0]->_run('execute');
}

sub _run {
    my $self = shift;
    my $mode = shift;

    for my $interval ($self->intervals) {
        $self->debug_message('Attempting to %s for interval %s.', $mode, $interval);
        my $cmd = Genome::Model::SingleSampleGenotype::Command::HaplotypeCaller->create(
            build => $self->build,
            intervals => [$interval],
        );

        unless($cmd->$mode) {
            $self->fatal_message('Could not %s for interval %s.', $mode, $interval);
        }

        $self->debug_message('Successfully %s for interval %s.', $mode, $interval);
    }

    return 1;
}

sub intervals {
    my $self = shift;
    my $build = $self->build;
    my $bucket = $self->bucket;

    unless ($bucket) {
        $self->fatal_message('No bucket supplied to bucket iterator');
    }

    return @{ $build->buckets_result->bucket($bucket) };
}

1;
