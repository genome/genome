package Genome::Config::Command::Validate;

use strict;
use warnings;

use Genome qw();
use Genome::Config qw();

class Genome::Config::Command::Validate {
    is => 'Command::V2',
};

sub execute {
    my $self = shift;

    my %env;
    my @all_specs = Genome::Config::all_specs();
    for my $spec (@all_specs) {
        my $value = Genome::Config::get($spec->key);

        if ($spec->has_env) {
            push @{ $env{$spec->env} }, $spec;
        }

        if ($spec->has_default_value) {
            my @errors = $spec->validate($spec->default_value);
            if (@errors) {
                my $msg = $spec->validation_error(@errors);
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
