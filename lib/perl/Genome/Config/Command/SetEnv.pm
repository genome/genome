package Genome::Config::Command::SetEnv;

use strict;
use warnings;

use Genome qw();
use Genome::Carp qw(dief);
use Genome::Config qw();

class Genome::Config::Command::SetEnv {
    doc => q(help set configuration variables via shell),
    is => 'Command::V2',
    has => [
        key => {
            is => 'Text',
            shell_args_position => 1,
            doc => 'The key for the configuration variable you wish to override.',
        },
        value => {
            is => 'Text',
            shell_args_position => 2,
            doc => 'The value you wish to use.',
        },
        format => {
            is => 'Text',
            is_optional => 1,
            default_value => 'bash',
            valid_values => [qw(bash tcsh)],
            doc => 'Specifies the output format.',
        },
    ],
};

sub help_detail {
    return <<EOS
Returns output that you can eval to override a configuration variable via the
shell environment.  For example, in a Bash script you might do,

    #!/bin/bash
    ...
    eval "\$(genome config set-env test_url http://localhost:8080/)"
    ...

EOS
}

sub execute {
    my $self = shift;

    if ($self->format eq 'bash') {
        print_bash($self->key, $self->value);
    }
    elsif ($self->format eq 'tcsh') {
        print_tcsh($self->key, $self->value);
    }
    else {
        dief 'unimplemented format: %s', $self->format;
    }

    return 1;
}

sub print_bash {
    my ($key, $value) = @_;
    my $spec = Genome::Config::spec($key);
    printf qq(export %s="%s"\n), $spec->env, $value;
}

sub print_tcsh {
    my ($key, $value) = @_;
    my $spec = Genome::Config::spec($key);
    printf qq(setenv %s "%s"\n), $spec->env, $value;
}

1;

