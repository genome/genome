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
            doc => 'Genome::Process',
        },
        destination => {
            is => 'Path',
            shell_args_position => 2,
            doc => 'The desired location of the created symlink(s)',
        }
    ],
};

sub help_synopsis {
    my $self = shift;
    my $result .= <<EOP;
    Creates symlinks to the results of a process.
EOP
    return $result;
}

sub help_detail {
    my $self = shift;
    my $result .= <<EOP;
The <destination> option should be an existing directory, or one that could be
created.

    genome process symlink-results <process_id> <destination>
EOP
    return $result;
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
