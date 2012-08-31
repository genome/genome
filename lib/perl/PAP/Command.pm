package PAP::Command;

use strict;
use warnings;

use PAP;
use Command;

use File::Temp;


class PAP::Command {
    is           => ['Command'],
    english_name => 'pap command',
};

sub command_name {

    my $self  = shift;
    my $class = ref($self) || $self;


    return 'pap' if $class eq __PACKAGE__;
    return $self->SUPER::command_name(@_);
    
}

sub help_brief {
    "modularized commands for testing Workflow"
}


1;
