package Genome::Command::Viewer;

use strict;
use warnings;

use Term::ReadKey 'GetTerminalSize';
use IO::Handle;

class Genome::Command::Viewer {
    doc => "Base class for 'view' commands.",
    is => 'Command::V2',
};


sub write_report {
    my ($self, $width, $handle) = @_;
    Carp::croak("Method write_report must be implemented in sub-classes.");
}

sub get_terminal_width {
    my $screen_width = 80;

    # this can fail in cases where no terminal is queriable
    eval {($screen_width) = GetTerminalSize();};

    return $screen_width;
}

sub execute {
    my ($self) = @_;


    my $handle = new IO::Handle;
    STDOUT->autoflush(1);
    $handle->fdopen(fileno(STDOUT), 'w');

    my $screen_width = $self->get_terminal_width();
    $self->write_report($screen_width, $handle);
    return 1;
}

sub get_report {
    my ($self, $width) = @_;

    my $handle = new IO::String;
    $self->write_report($width, $handle);
    my $report = ${$handle->string_ref};
    $handle->close();

    return $report;
}

1;
