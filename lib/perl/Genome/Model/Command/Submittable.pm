package Genome::Model::Command::Submittable;

use strict;
use warnings;

use Genome;

use Genome::Sys::LSF::bsub qw();

role Genome::Model::Command::Submittable {
};

sub _submit_jobs {
    my $self = shift;
    my $anp = shift;
    my $cmd = shift;

    my @vol = $anp->possible_volumes;

    my $anp_guard = $anp->set_env;
    my $volume_guard = Genome::Config::set_env('docker_volumes', join(' ', map { "$_:$_" } @vol));
    my $image_guard = Genome::Config::set_env('lsb_sub_additional', sprintf('docker0(%s)', $ENV{LSB_DOCKER_IMAGE})); #use the current image regardless of the AnP config
    unless($ENV{LSF_DOCKER_NETWORK} eq 'host') {
        $self->fatal_message('Parent container must have LSF_DOCKER_NETWORK=host set.');
    }

    local $ENV{UR_NO_REQUIRE_USER_VERIFY} = 1;

    Genome::Sys::LSF::bsub::bsub(
        cmd => $cmd,
        user_group => Genome::Config::get('lsf_user_group'),
        interactive => 1,
        queue => Genome::Config::get('lsf_queue_interactive'),
    );

    return 1;
}

1;
