package Genome::Model::Command::Admin::PurgeDiskAllocationsForAbandonedBuilds;

use strict;
use warnings;

use Genome;

class Genome::Model::Command::Admin::PurgeDiskAllocationsForAbandonedBuilds {
    is => 'Command::V2',
    roles => ['Genome::Model::Command::Submittable'],
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

        for my $build_class (keys %{$builds_by_anp{$anp_id}}) {
            my @build_ids = keys %{ $builds_by_anp{$anp_id}{$build_class} };
            my $connector = (scalar(@build_ids) == 1)? '=' : ':';

            $self->_submit_jobs(
                $anp,
                [
                    qw(genome disk allocation purge --reason),
                    'purging leftover allocation for abandoned build',
                    'status=active,owner_class_name=' . $build_class . ',owner_id' . $connector . join("/", @build_ids),
                ],
            );
        }
    }

    return 1;
}

1;
