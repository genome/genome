package Genome::Model::Build::Command::SymlinkResults;

use strict;
use warnings;

use Genome;

class Genome::Model::Build::Command::SymlinkResults {
    is => ['Command::V2'],
    has_input => [
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
    has_transient_output => [
        output_directory => {
            is => 'DirectoryPath',
            doc => 'The directory of symlinks that was created',
        },
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

    my $dir = $build->symlink_results($destination);
    $self->output_directory($dir);

    return 1;
}

1;
