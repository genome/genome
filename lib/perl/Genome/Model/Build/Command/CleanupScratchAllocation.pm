package Genome::Model::Build::Command::CleanupScratchAllocation;

use strict;
use warnings;

use Genome;
use Try::Tiny qw(try catch);

class Genome::Model::Build::Command::CleanupScratchAllocation {
    is  => 'Genome::Model::Build::Command::Base',
};

sub help_brief {
    return "Cleanup the scratch allocation for a build";
}

sub help_detail {
    return "This command cleans up the scratch allocation ('tmp') directory for a build while otherwise leaving its data directory intact.  This can be useful if the logs for a failed build need to preserved for further inspection, but the working space needs to be freed."
}

sub execute {
    my $self = shift;

    my @builds = $self->builds;
    for my $build (@builds) {
         $build->cleanup_scratch_directory;
    }

    return 1;
}

1;
