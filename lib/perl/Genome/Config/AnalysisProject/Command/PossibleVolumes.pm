package Genome::Config::AnalysisProject::Command::PossibleVolumes;

use strict;
use warnings;

use feature qw(say);

use Genome;

class Genome::Config::AnalysisProject::Command::PossibleVolumes {
    is => 'Genome::Config::AnalysisProject::Command::Base',
    has_optional => [
        use_docker_format => {
            is => 'Boolean',
            default => 0,
            doc => 'If set, will output in a format suitable for the "-v" option to `docker run`.',
        },
    ],
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

    my @possible_volumes = $self->analysis_project->possible_volumes;

    if ($self->use_docker_format) {
        say join(" ", map { "$_:$_" } @possible_volumes);
    } else {
        say $_ for @possible_volumes;
    }

    return 1;
}

1;
