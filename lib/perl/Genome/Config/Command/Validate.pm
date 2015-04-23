package Genome::Config::Command::Validate;

use strict;
use warnings;

use Genome qw();
use Genome::Config qw();

class Genome::Config::Command::Validate {
    doc => 'validate configuration key-value pairs',
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
    'Validate configuration key-value pairs.'
}

sub __default_keys__ {
    return [ Genome::Config::all_keys() ];
}

sub execute {
    my $self = shift;

    my %env;
    for my $key ($self->keys) {
        my $value = Genome::Config::get($key);

        my $spec = Genome::Config::spec($key);

        if ($spec->has_env) {
            push @{ $env{$spec->env} }, $spec;
        }

        if ($spec->has_default_value) {
            my @errors = $spec->validate($spec->default_value);
            if (@errors) {
                my $msg = Genome::Config::validation_error($spec, @errors);
                croakf($msg);
            }
        }
    }

    for my $env_key (keys %env) {
        my @specs = @{ $env{$env_key} };
        if (@specs > 1) {
            my $spec_keys = join(',', map { $_->key } @specs);
            croakf('%s is used as the env value for multiple keys: %s', $env_key, $spec_keys);
        }
    }

    return 1;
}

1;
