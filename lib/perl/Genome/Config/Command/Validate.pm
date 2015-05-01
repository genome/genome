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

    my $return = 1;
    my %env;
    for my $key ($self->keys) {
        my $spec = Genome::Config::spec($key);

        my $error = Genome::Config::validate($key);
        if (defined $error) {
            my $msg = $spec->validation_error($error);
            printf("%s\n", $msg);
            $return = 0;
        }

        if ($spec->has_env) {
            push @{ $env{$spec->env} }, $spec;
        }

        if ($spec->has_default_value) {
            my $error = $spec->validate($spec->default_value);
            if (defined $error) {
                printf("%s has an invalid default_value\n", $spec->key);
                $return = 0;
            }
        }
    }

    for my $env_key (keys %env) {
        my @specs = @{ $env{$env_key} };
        if (@specs > 1) {
            my $spec_keys = join(',', map { $_->key } @specs);
            printf("%s is used as the env value for multiple keys: %s\n", $env_key, $spec_keys);
            $return = 0;
        }
    }

    return $return;
}

1;
