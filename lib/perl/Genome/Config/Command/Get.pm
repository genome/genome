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
        format => {
            is => 'Text',
            is_optional => 1,
            default_value => 'default',
            valid_values => [qw(bash default tcsh)],
            doc => 'Specifies the output format.',
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
        $self->print($key);
    }
    return 1;
}

sub print {
    my ($self, $key) = @_;

    if ($self->format eq 'bash') {
        return print_bash($key);
    }

    if ($self->format eq 'tcsh') {
        return print_tcsh($key);
    }

    return print_default($key);
}

sub print_bash {
    my $key = shift;
    my $spec = Genome::Config::spec($key);
    my $value = Genome::Config::get($key);
    printf qq(%s="%s"\n), $spec->env, $value;
}

sub print_default {
    my $key = shift;
    my $value = Genome::Config::get($key);
    printf "%s = '%s'\n", $key, $value;
}

sub print_tcsh {
    my $key = shift;
    my $spec = Genome::Config::spec($key);
    my $value = Genome::Config::get($key);
    printf qq(setenv %s "%s"\n), $spec->env, $value;
}

1;
