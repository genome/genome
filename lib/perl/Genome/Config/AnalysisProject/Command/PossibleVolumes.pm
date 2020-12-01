package Genome::Config::AnalysisProject::Command::PossibleVolumes;

use strict;
use warnings;

use feature qw(say);

use Genome;

class Genome::Config::AnalysisProject::Command::PossibleVolumes {
    is => 'Genome::Config::AnalysisProject::Command::Base',
};

sub help_brief {
    return 'output the possible mount paths for the disk groups in the environment config';
}

sub help_detail {
    return <<EOS
This command is used to list all the possible volume mount paths that might be used for an analysis-project based on its environment configuration.  This forms the set of volumes that should be mounted inside any containers that need to do work for the analysis-project.
EOS
}

sub valid_statuses {
    return ('Pending', 'In Progress', 'Hold');
}

sub execute {
    my $self = shift;

    my $guard = $self->analysis_project->set_env;

    my @group_names = map { Genome::Config::get($_) } (qw(disk_group_models disk_group_alignments disk_group_scratch));
    my @volumes = Genome::Disk::Volume->get(disk_group_names => \@group_names);

    for my $v (sort { $a->mount_path cmp $b->mount_path } @volumes) {
        say $v->mount_path;
    }

    return 1;
}

1;
