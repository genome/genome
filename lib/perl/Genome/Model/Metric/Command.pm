package Genome::Model::Metric::Command;

use strict;
use warnings;

use Genome;

class Genome::Model::Metric::Command {
    is => 'Command',
    doc => "work with per-build model metrics",
};

sub command_name {
    my $class = ref($_[0]) || $_[0];
    return $class->SUPER::command_name unless $class eq __PACKAGE__;
    return 'genome model metric';
}

sub command_name_brief {
    my $class = ref($_[0]) || $_[0];
    return $class->SUPER::command_name_brief unless $class eq __PACKAGE__;
    return 'metric';
}

1;

