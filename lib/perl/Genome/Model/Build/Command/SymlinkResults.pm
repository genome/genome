package Genome::Model::Build::Command::SymlinkResults;

use strict;
use warnings;

use Genome;

class Genome::Model::Build::Command::SymlinkResults {
    is => ['Command::V2'],
    has => [
        build => {
            is => 'Genome::Model::Build',
            shell_args_position => 1,
            doc => "The builds that has results you'd like to create symlinks to",
        },
        destination => {
            is => 'DirectoryPath',
            shell_args_position => 2,
            doc => 'The desired location under which to make a directory of symlink(s)',
        }
    ],
    doc => 'Creates symlinks to the results of a build.',
};

sub help_detail {
    my $self = shift;
    return <<EOP;
Creates symlinks to the results of a build.  The <destination> option should
be an existing directory.
EOP
}

sub execute {
    my $self = shift;
    my $build = $self->build;
    my $destination = $self->destination;

    Genome::Sys->validate_directory_for_read_write_access($destination);

    unless($build->can('symlink_results')) {
        $self->fatal_message('Symlinking results is unavailable for builds of type "%s"', $build->type_name);
    }

    my $dir = $build->symlink_results($destination);
    $self->status_message('Results for build %s symlinked to %s', $build->__display_name__, $dir);

    return 1;
}

1;
