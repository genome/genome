package GAP::Command;

use strict;
use warnings;

use GAP;
use Command;

class GAP::Command {
    is => ['Command'],
    english_name => 'gap command',
};

sub command_name {
    my $self = shift;
    my $class = ref($self) || $self;
    return 'gap' if $class eq __PACKAGE__;
    return $self->SUPER::command_name(@_);
}

sub help_brief {
    "modularized commands for testing Workflow"
}

1;
