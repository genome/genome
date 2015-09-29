package Genome::Ptero::TestCommand::OptionalInputCommand;

use strict;
use warnings;
use Genome;

class Genome::Ptero::TestCommand::OptionalInputCommand {
    is => "Command::V2",
    has_input => [
        required_input => {
            is => "Text",
        },
        optional_input => {
            is => "Text",
            is_optional => 1,
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

    if (defined($self->optional_input)) {
        $self->output($self->optional_input)
    } else {
        $self->output($self->required_input)
    }
    return 1;
};

1;
