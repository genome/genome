package Genome::Model::Command::Admin::PurgeDiskAllocationsForAbandonedBuilds;

use strict;
use warnings;

use Genome;

class Genome::Model::Command::Admin::PurgeDiskAllocationsForAbandonedBuilds {
    is => 'Command::V2',
    doc => 'Purge remaining allocations for abandoned builds.',
    has => [
        submit_purge_jobs => {
            is => 'Boolean',
            default => 0,
            doc => 'if set, will use `genome analysis-project possible-volumes` to submit jobs for purging',
        },
    ],
    has_optional_transient => [
        builds => {
            is => 'Genome::Model::Build',
            require_user_verify => 1,
            is_many => 1,
        },
    ],
};

sub help_detail {
    'Purge active disk allocaitons remaining for abandoned builds.';
};

sub execute {
    my $self = shift;

    my @builds = $self->resolve_param_value_from_cmdline_text(
        {
            name => 'builds',
            class => 'Genome::Model::Build',
            value => ['status=Abandoned,run_by=' . Genome::Sys->username . ',disk_allocations.status=active'],
        }
    );

    if ($self->submit_purge_jobs) {
        $self->_submit_purge_jobs(@builds);
    } else {
        for my $build (@builds) {
            map $_->purge('purging leftover allocation for abandoned build'), $build->disk_allocations;
        }
    }
}

sub _submit_purge_jobs {
    my $self = shift;
    my @builds = @_;

    my %builds_by_anp;
    for my $build (@builds) {
        my $anp = $build->model->analysis_project;
        $builds_by_anp{$anp->id}{$build->class}{$build->id} = 1;
    }

    for my $anp_id (keys %builds_by_anp) {
        my $anp = Genome::Config::AnalysisProject->get($anp_id);

        my @vol = $anp->possible_volumes;

        my $anp_guard = $anp->set_env;
        my $volume_guard = Genome::Config::set_env('docker_volumes', join(' ', map { "$_:$_" } @vol));
        my $image_guard = Genome::Config::set_env('lsb_sub_additional', sprintf('docker0(%s)', $ENV{LSB_DOCKER_IMAGE})); #use the current image regardless of the AnP config
        unless($ENV{LSF_DOCKER_NETWORK} eq 'host') {
            $self->fatal_message('Parent container must have LSF_DOCKER_NETWORK=host set.');
        }

        for my $build_class (keys %{$builds_by_anp{$anp_id}}) {

            local $ENV{UR_NO_REQUIRE_USER_VERIFY} = 1; 

            my @build_ids = keys %{ $builds_by_anp{$anp_id}{$build_class} };
            my $connector = (scalar(@build_ids) == 1)? '=' : ':';

            Genome::Sys::LSF::bsub::bsub(
                cmd => [
                    qw(genome disk allocation purge --reason),
                    'purging leftover allocation for abandoned build',
                    'status=active,owner_class_name=' . $build_class . ',owner_id' . $connector . join("/", @build_ids),
                ],
                user_group => Genome::Config::get('lsf_user_group'),
                interactive => 1,
            );
        }
    }

    return 1;
}

1;
