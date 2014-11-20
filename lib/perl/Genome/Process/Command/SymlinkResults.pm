package Genome::Process::Command::SymlinkResults;

use strict;
use warnings FATAL => 'all';
use Genome;

class Genome::Process::Command::SymlinkResults {
    is => ['Command::V2'],
    has => [
        process => {
            is => 'Genome::Process',
            shell_args_position => 1,
            doc => "The Genome::Process that has results you'd like to create symlinks to",
        },
        destination => {
            is => 'Path',
            shell_args_position => 2,
            doc => 'The desired location of the created symlink(s)',
        }
    ],
    doc => 'Creates symlinks to the results of a process.',
};

sub help_detail {
    my $self = shift;
    return <<EOP;
Creates symlinks to the results of a process.  The <destination> option should
be an existing directory, or one that could be created.
EOP
}

sub shortcut {
    my $self = shift;
    return $self->execute();
}

sub execute {
    my $self = shift;
    my $rv = $self->process->symlink_results($self->destination);
    if ($rv) {
        $self->status_message("Successfully symlinked results for process (%s) to (%s)",
            $self->process->id, $self->destination);
    }
    return $rv;
}


1;
