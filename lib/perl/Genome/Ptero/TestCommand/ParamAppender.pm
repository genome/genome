package Genome::Ptero::TestCommand::ParamAppender;

use strict;
use warnings;
use Genome;

class Genome::Ptero::TestCommand::ParamAppender {
    is => "Command::V2",
    has_input => [
        prefix => {
            is => "Text",
        },
        suffix => {
            is => "Text",
        },
    ],
    has_output => [
        output => {
            is => "Text",
        },
    ],
};

sub shortcut {
    my $self = shift;
    return $self->execute(@_);
}

sub execute {
    my $self = shift;
    my $p = $self->prefix;
    my $s = $self->suffix;
    my $a = "$p$s";
    print "$p + $s = $a!!\n";
    $self->output($a);
    return 1;
};

1;
