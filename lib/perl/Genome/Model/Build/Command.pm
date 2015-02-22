package Genome::Model::Build::Command;

use strict;
use warnings;

use Genome;

class Genome::Model::Build::Command {
    is => 'Command',
    doc => "work with model builds",
    has => [
        builds => {
            is => 'Genome::Model::Build',
            is_many => 1,
            shell_args_position => 1,
        },

    ],
};

#< Name >#
sub command_name {
    my $class = ref($_[0]) || $_[0];
    return $class->SUPER::command_name unless $class eq __PACKAGE__;
    return 'genome model build';
}

sub command_name_brief {
    my $class = ref($_[0]) || $_[0];
    return $class->SUPER::command_name_brief unless $class eq __PACKAGE__;
    return 'build';
}
#<>#

1;

