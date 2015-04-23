package Genome::Config::Command::Get;

use strict;
use warnings;

use Genome qw();
use Genome::Config qw();

class Genome::Config::Command::Get {
    doc => 'list configuration key-value pairs',
    is => 'Command::V2',
    has => [
        keys => {
            is => 'Text',
            is_many => 1,
            is_optional => 1,
            shell_args_position => 1,
            calculated_default => 1,
            doc => 'Limit to one or more keys; otherwise all keys.',
        },
    ],
};

sub help_detail {
    'List configuration key-value pairs.'
}

sub __default_keys__ {
    return [ Genome::Config::all_keys() ];
}

sub execute {
    my $self = shift;
    for my $key ($self->keys) {
        my $value = Genome::Config::get($key);
        printf "%s = '%s'\n", $key, $value;
    }
    return 1;
};

1;
