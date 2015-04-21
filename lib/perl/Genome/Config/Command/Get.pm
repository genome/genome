package Genome::Config::Command::Get;

use strict;
use warnings;

use Genome qw();
use Genome::Config qw();

class Genome::Config::Command::Get {
    is => 'Command::V2',
    has => [
        key => {
            is => 'Text',
            shell_args_position => 1,
        },
    ],
};

sub execute {
    my $self = shift;
    my $value = Genome::Config::get($self->key);
    print $value, "\n";
    return 1;
};

1;
