package Genome::Config::Command::Get;

use strict;
use warnings;

use Genome;
use Genome::Config qw();

class Genome::Config::Command::Get {
    doc => 'get the value for a specific configuration variable',
    is => 'Command::V2',
    has => [
        key => {
            is => 'Text',
            shell_args_position => 1,
            doc => "The key of the configuration variable."
        },
    ],
};

sub help_detail {
    return <<EOS
Get the value for a specific configuration variable (by key).  For example, in
a Bash script you might do,

    #!/bin/bash
    ...
    TEST_URL="\$(genome config get test_url)"
    ...

EOS
}

sub execute {
    my $self = shift;
    my $value = Genome::Config::get($self->key);
    printf "%s\n", $value;
    return 1;
}

1;
