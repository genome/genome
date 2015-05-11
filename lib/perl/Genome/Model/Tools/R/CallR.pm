package Genome::Model::Tools::R::CallR;

use warnings;
use strict;
use Genome;
use Cwd;
use Statistics::R;
require Genome::Sys;

class Genome::Model::Tools::R::CallR {
    is => 'Command',
    has => [
    command => {
        type => 'String',
        is_optional => 0,
        is_input => 1,
        doc => 'R command to be run in temp directory using this script',
    },
    library => {
        type => 'String',
        is_optional => 0,
        is_input => 1,
        doc => 'R library containing functions for specific application. Something like "cmds_lib.R". See the GMT/R directory for library files labeled like this: CallR.pm.something_lib.R, where "something_lib.R" is the input to this option.',
    },
    r_session_info => {
        type => 'Boolean',
        is_optional => 1,
        default => 0,
        doc => 'flag to get R session info.',
    },
    ]
};

sub help_brief {
    "Wrapper to call functions in R library."
}
sub help_detail {
    "Wrapper to call functions in R library."
}

sub execute {
    my $self = shift;

    my $command = $self->command;
    my $r_library = __FILE__ . "." . $self->library;
    
    my $tempdir = Genome::Sys->create_temp_directory();
    my $cwd = cwd();
    my $R = Statistics::R->new(tmp_dir => $tempdir);
    $R->startR();
    if($self->r_session_info) {
        $self->print_session_info($R);
    }
    $R->send("source('$r_library');");
    $R->send($command);
    $R->stopR();
    chdir $cwd;
    return 1;
}

sub print_session_info {
    my $self = shift;
    my $R = shift;
    $R->send("print(sessionInfo())");
    $self->status_message("SessionInfo\n" . $R->read(20) . "\n");
    $R->send("print(.libPaths())");
    $self->status_message("LibPaths\n" . $R->read(20) . "\n");
    $R->send("print(capabilities())");
    $self->status_message("Capabilities\n" . $R->read(20) . "\n");
}

1;
