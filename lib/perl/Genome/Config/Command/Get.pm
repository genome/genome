package Genome::Config::Command::Get;

use strict;
use warnings;

use Genome qw();
use Genome::Config qw();

class Genome::Config::Command::Get {
    is => 'Command::V2',
    has => [
        keys => {
            is => 'Text',
            is_many => 1,
            is_optional => 1,
            shell_args_position => 1,
        },
    ],
};

sub execute {
    my $self = shift;
    my @keys = $self->keys || Genome::Config::all_keys();
    for my $key (@keys) {
        my $value = Genome::Config::get($key);
        printf "%s = '%s'\n", $key, $value;
    }
    return 1;
};

1;
